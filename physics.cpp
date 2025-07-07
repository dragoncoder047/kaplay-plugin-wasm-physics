#include <emscripten.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <optional>
#include <limits>

#define EXPORT EMSCRIPTEN_KEEPALIVE extern "C"

typedef double f64;
#define EPSILON (std::numeric_limits<f64>::epsilon())

// MARK: imports

/** handles the collisions, when found */
EM_IMPORT(handleCollision) void handleCollision(int idA, int idB, f64 normX, f64 normY, f64 distance);

/** returns true if the objects shouldn't collide with each other. */
EM_IMPORT(checkCollisionIgnore) bool checkCollisionIgnore(int idA, int idB);

/** returns true if the object isn't paused */
EM_IMPORT(isObjActive) bool isObjActive(int id);

/* tells the JS code to send the transformation data for the object with this id */
EM_IMPORT(pullTransform) void pullTransform(int id);

// MARK: lists n' stuff

typedef void (*iter_body)(void *arg, void *item);
typedef f64 (*comparator)(void *left, void *right);

class List {
    public:
    class Node {
        public:
        void *item;
        Node *next;
        Node(void *item, Node *next): item(item), next(next) {}
    };

    Node *head;
    int length;

    List(): head(NULL), length(0) {}
    ~List() { clear(); }
    void clear() {
        while (head != NULL) {
            Node *next = head->next;
            delete head;
            head = next;
        }
        length = 0;
    }
    void push(void *item) { head = new Node(item, head); length++; }
    void iterate(void *arg, iter_body body) {
        for (Node *node = head; node != NULL; node = node->next)
            body(arg, node->item);
    }
    void *remove(void *target, comparator found) {
        for (Node **node = &head; *node != NULL; node = &(*node)->next) {
            if (found(target, (*node)->item)) {
                Node *removed = *node;
                void *item = removed->item;
                *node = (*node)->next;
                delete removed;
                length--;
                return item;
            }
        }
        return NULL;
    }
    void *find(void *target, comparator found) {
        for (Node *node = head; node != NULL; node = node->next) {
            if (found(target, node->item)) {
                return node->item;
            }
        }
        return NULL;
    }
    void reverse() {
        Node *previous = NULL;
        Node *current = head;
        Node *next = NULL;
        while (current != NULL) {
            next = current->next;
            current->next = previous;
            previous = current;
            current = next;
        }
        head = previous;
    }
    void insertionSort(comparator compare) {
        Node *sorted = NULL;
        while (head != NULL) {
            Node *current = head;
            head = head->next;
            Node **p = &sorted;
            while (*p != NULL && compare(current->item, (*p)->item) < 0.)
                p = &((*p)->next);
            current->next = *p;
            *p = current;
        }
        head = sorted;
    }
};

// MARK: mathy math math

class vec2 {
    public:
    f64 x, y;
    vec2(): x(0.), y(0.) {}
    vec2(f64 x, f64 y): x(x), y(y) {}
    f64 length() const { return sqrt(x*x + y*y); }
    inline vec2 operator+(vec2 other) const { return vec2(x + other.x, y + other.y); }
    inline vec2 operator-(vec2 other) const { return vec2(x - other.x, y - other.y); }
    inline vec2 operator-() const { return vec2(-x, -y); }
    inline vec2 operator*(f64 other) const { return vec2(x * other, y * other); }
    inline vec2 operator*(vec2 other) const { return vec2(x * other.x, y * other.y); }
    inline vec2 operator/(f64 other) const { return vec2(x / other, y / other); }
    inline vec2 unit() const { f64 l = length(); return vec2(x/l, y/l); }
    inline f64 dot(vec2 other) const { return x * other.x + y * other.y; }
    inline vec2 rotate(f64 angle) const { f64 s = sin(angle), c = cos(angle); return vec2(x*c - y*s, x*s + y*c); }
};

class mat23 {
    public:
    f64 a, b, c, d, e, f;
    mat23(): a(1.), b(0.), c(0.), d(1.), e(0.), f(0.) {}
    mat23(f64 a, f64 b, f64 c, f64 d, f64 e, f64 f):
        a(a), b(b), c(c), d(d), e(e), f(f) {}
    inline mat23 operator+(vec2 translate) const { return mat23(a, b, c, d, e + vec2(a, c).dot(translate), f + vec2(b, d).dot(translate)); }
    inline mat23 operator*(vec2 scale) const { return mat23(a * scale.x, b * scale.x, c * scale.y, d * scale.y, e, f); }
    inline vec2 transformVector(vec2 v) const { return vec2(vec2(a, c).dot(v), vec2(b, d).dot(v)); }
    inline vec2 transformPoint(vec2 v) const { return vec2(e, f) + transformVector(v); }
    inline f64 getRotation() const { return atan2(-c, a); }
    inline vec2 getScale() const { return vec2(vec2(a, c).length(), vec2(b, d).length()); }
};

typedef struct { f64 e1, e2; } eigenvalues_t;
typedef struct { vec2 v1, v2; } eigenvectors_t;

class mat2 {
    public:
    f64 a, b, c, d;
    mat2(): a(1.), b(0.), c(0.), d(1.) {}
    mat2(f64 a, f64 b, f64 c, f64 d): a(a), b(b), c(c), d(d) {}
    f64 det() const { return a * d - b * c; }
    mat2 inverse() const { f64 dt = det(); return mat2(d/dt, -b/dt, -c/dt, a/dt); }
    mat2 transpose() const { return mat2(a, c, b, d); }
    f64 trace() const { return a + d; }
    mat2 operator*(mat2 o) const { return mat2(a * o.a, b * o.b, c * o.c, d * o.d); }
    vec2 transform(vec2 p) const { return vec2(vec2(a, b).dot(p), vec2(c, d).dot(p)); }
    eigenvalues_t eigenvalues() const {
        f64 m = trace() / 2., d = det();
        f64 t = sqrt(m * m - d);
        return {m + t, m - t};
    }
    eigenvectors_t eigenvectors(eigenvalues_t ev) const {
        if (c != 0) return {vec2(ev.e1 - d, c), vec2(ev.e2 - d, c)};
        if (b != 0) return {vec2(b, ev.e1 - a), vec2(b, ev.e2 - a)};
        if (abs(transform(vec2(1., 0.)).x - ev.e1) < EPSILON)
            return {vec2(1., 0.), vec2(0., 1.)};
        return {vec2(0., 1.), vec2(1., 0.)};
    }
};

class mat3 {
    public:
    f64 m11, m12, m13, m21, m22, m23, m31, m32, m33;
    mat3(f64 m11, f64 m12, f64 m13, f64 m21, f64 m22, f64 m23, f64 m31, f64 m32, f64 m33):
        m11(m11), m12(m12), m13(m13), m21(m21), m22(m22), m23(m23), m31(m31), m32(m32), m33(m33) {}
    mat3(mat2 f): m11(f.a), m12(f.b), m13(0.), m21(f.c), m22(f.d), m23(0.), m31(0.), m32(0.), m33(1.) {}
    inline mat2 toMat2() const { return mat2(m11, m12, m21, m22); }
    inline mat3 rotate(f64 angle) const {
        f64 c = cos(angle), s = sin(angle);
        return mat3(c * m11 + s * m21, c * m12 + s * m22, m13,
                    c * m21 - s * m11, c * m22 - s * m12, m23,
                    m31,               m32,               m33);
    }
    inline mat3 scale(vec2 s) const {
        return mat3(m11 * s.x, m12 * s.x, m13,
                    m21 * s.y, m22 * s.y, m23,
                    m31,       m32,       m33);
    }
};

// MARK: GJK algorithm

typedef struct { f64 left, right; } EdgeBounds;

enum ColliderType {
    NONE,
    ELLIPSE,
    POLYGON
};
EXPORT ColliderType CT_ELLIPSE = ColliderType::ELLIPSE;
EXPORT ColliderType CT_POLYGON = ColliderType::POLYGON;

class Collider {
    public:
    vec2 center;
    vec2 anchor;
    ColliderType type;
    Collider(): type(NONE) {}
    virtual ~Collider() {} // allow deletion by any pointer
    Collider(vec2 center, ColliderType type):
        center(center), anchor(vec2()), type(type) {}
    virtual vec2 support(vec2 direction) const = 0;
    virtual EdgeBounds bounds() const = 0;
    virtual void transform(mat23 t, Collider **out) = 0;
};

class Ellipse: public Collider {
    public:
    f64 angle;
    vec2 radii;
    Ellipse(vec2 center, vec2 radii, f64 angle):
        Collider(center, ELLIPSE), angle(angle), radii(radii) {}
    vec2 support(vec2 direction) const {
        return (center
            + (angle == 0
                // Axis aligned
                ? (direction.unit() * radii)
                // Rotated
                : (direction.rotate(-angle).unit() * radii).rotate(angle)));
    }
    EdgeBounds bounds() const {
        if (angle == 0) {
            return {center.x - radii.x, center.x + radii.x};
        }
        f64 halfwidth = (radii * vec2(1, 0).rotate(angle)).length();
        return {center.x - halfwidth, center.x + halfwidth};
    }
    void transform(mat23 t, Collider **out) {
        if (out == NULL) abort();
        if (*out == NULL || (*out)->type != ELLIPSE) {
            delete *out;
            *out = new Ellipse(center, radii, angle);
        }
        // TODO: should this also take anchor into account
        Ellipse *res = (Ellipse *)*out;
        res->center = t.transformPoint(center);
        if (angle == 0 && t.getRotation() == 0) {
            res->radii = vec2(t.a, t.d) * radii;
            res->angle = 0;
        } else {
            res->fromMat2(mat3(toMat2()).scale(t.getScale()).rotate(t.getRotation()).toMat2());
        }
    }
    inline mat2 toMat2() const {
        f64 c = cos(angle), s = sin(angle);
        return mat2(c * radii.x, -s * radii.y, s * radii.x, c * radii.y);
    }
    inline void fromMat2(mat2 tr) {
        mat2 inv = tr.inverse();
        mat2 M = inv.transpose() * inv;
        eigenvalues_t eva = M.eigenvalues();
        eigenvectors_t evc = M.eigenvectors(eva);

        f64 a = 1/sqrt(eva.e1), b = 1/sqrt(eva.e2);

        // Make sure we use the semi-major axis for the rotation
        if (a > b) {
            radii = vec2(a, b);
            angle = atan2(-evc.v1.y, evc.v1.x);
        }
        else {
            radii = vec2(b, a);
            angle = atan2(-evc.v2.y, evc.v2.x);
        }
    }
};

class Polygon: public Collider {
    public:
    vec2 *vertices = NULL;
    bool isRect = false;
    int length = 0, capacity = 0;
    Polygon(vec2 center): Collider(center, POLYGON) { addVertex(center); }
    ~Polygon() { delete[] vertices; }
    void ensureCapacity(int length) {
        if (length <= capacity) return;
        while (length > capacity) capacity <<= 1;
        vec2 *newVertices = new vec2[capacity];
        memcpy(newVertices, vertices, length * sizeof(vec2));
        delete[] vertices;
        vertices = newVertices;
    }
    void addVertex(vec2 v) {
        ensureCapacity(length + 1);
        vertices[length++] = v;
    }
    vec2 support(vec2 direction) const {
        vec2 maxPoint = vec2(NAN, NAN);
        f64 maxDistance = -INFINITY;
        for (size_t i = 0; i < length; i++) {
            vec2 vertex = vertices[i];
            f64 dist = vertex.dot(direction);
            if (dist > maxDistance) {
                maxDistance = dist;
                maxPoint = vertex;
            }
        }
        return maxPoint;
    }
    EdgeBounds bounds() const {
        f64 min = INFINITY, max = -INFINITY;
        for (size_t i = 0; i < length; i++) {
            vec2 vertex = vertices[i];
            if (vertex.x > max) max = vertex.x;
            if (vertex.x < min) min = vertex.x;
        }
        return {min, max};
    }
    void transform(mat23 t, Collider **out) {
        if (out == NULL) abort();
        if (*out == NULL || (*out)->type != POLYGON) {
            delete *out;
            *out = new Polygon(vec2(0, 0));
        }
        if (isRect) {
            if (length != 4) abort();
            vec2 dims = vertices[2] - vertices[0];
            t = t + ((anchor + vec2(1., 1.)) * dims * -.5);
        }
        Polygon *res = (Polygon *)*out;
        res->ensureCapacity(length);
        for (int i = 0; i < length; i++) {
            res->vertices[i] = t.transformPoint(vertices[i]);
        }
        res->center = res->vertices[0];
    }
};

inline vec2 calculateSupport(Collider *a, Collider *b, vec2 direction) {
    // Calculate the support vector. This is done by calculating the difference between
    // the furthest points found of the shapes along the given direction.
    return a->support(direction) - b->support(-direction);
}

inline vec2 tripleProduct(vec2 a, vec2 b, vec2 c) {
    // AxB = (0, 0, axb)
    // AxBxC = (-axb * c.y, axb * c.x, 0)
    f64 n = a.x * b.y - a.y * b.x;
    // This vector lies in the same plane as a and b and is perpendicular to c
    return vec2(-n * c.y, n * c.x);
}

#define MAX_TRIES 20

typedef enum {
    EVOLVE_DONE_NO_INTERSECTION,
    EVOLVE_DONE_FOUND_INTERSECTION,
    EVOLVE_IN_PROGRESS,
} EvolveResult;

typedef enum {
    CLOCKWISE,
    COUNTERCLOCKWISE,
} PolygonWinding;

class SimplexEdge {
    public:
    f64 distance;
    int index;
    vec2 normal;
    SimplexEdge(f64 distance, int index, vec2 normal): distance(distance), index(index), normal(normal) {}
};

class GJKResult {
    public:
    vec2 normal;
    f64 distance;
    GJKResult(vec2 normal, f64 distance): normal(normal), distance(distance) {}
};

class Simplex {
    public:

    vec2 vertices[3];
    int length = 0;

    Simplex() {}

    void reset() { length = 0; }

    bool addSupport(Collider *a, Collider *b, vec2 direction) {
        vec2 support = calculateSupport(a, b, direction);
        vertices[length++] = support;
        // Returns true if both vectors are in the same direction
        return direction.dot(support) >= 0;
    }

    EvolveResult evolve(Collider *colliderA, Collider *colliderB, vec2 *direction) {
        switch (length) {
            case 0: {
                // Zero points, set the direction the center of colliderA
                // towards the center of of colliderB
                *direction = colliderB->center - colliderA->center;
                break;
            }
            case 1: {
                // Reverse the direction, to make a line
                *direction = -(*direction);
                break;
            }
            case 2: {
                // We now have a line ab. Take the vector ab and the vector a origin
                vec2 ab = vertices[1] - vertices[0];
                vec2 a0 = -vertices[0];

                // Get the vector perpendicular to ab and a0
                // Then get the vector perpendicular to the result and ab
                // This is our new direction to form a triangle
                *direction = tripleProduct(ab, a0, ab);
                break;
            }
            case 3:
                {
                    // We have a triangle, and need to check if it contains the origin
                    vec2 c0 = -vertices[2];
                    vec2 bc = vertices[1] - vertices[2];
                    vec2 ca = vertices[0] - vertices[2];

                    vec2 bcNorm = tripleProduct(ca, bc, bc);
                    vec2 caNorm = tripleProduct(bc, ca, ca);

                    if (bcNorm.dot(c0) > 0) {
                        // The origin does not lie within the triangle
                        // Remove the first point and look in the direction of bcNorm
                        vertices[0] = vertices[1];
                        vertices[1] = vertices[2];
                        length--;
                        *direction = bcNorm;
                    }
                    else if (caNorm.dot(c0) > 0) {
                        // The origin does not lie within the triangle
                        // Remove the second point and look in the direction of caNorm
                        vertices[1] = vertices[2];
                        length--;
                        *direction = caNorm;
                    }
                    else {
                        // The origin lies within the triangle
                        return EVOLVE_DONE_FOUND_INTERSECTION;
                    }
                }
                break;
            default:
                // too many vertices
                abort();
        }

        // Try to add a new support point to the simplex
        // If successful, continue evolving
        return addSupport(colliderA, colliderB, *direction)
            ? EVOLVE_IN_PROGRESS
            : EVOLVE_DONE_NO_INTERSECTION;
    }


    /** Returns the edge closest to the origin */
    SimplexEdge findClosestEdge(PolygonWinding winding) {
        f64 minDistance = INFINITY;
        vec2 minNormal, line, norm;
        int minIndex = 0;
        for (int i = 0; i < length; i++) {
            int j = (i + 1) % length;
            line = vertices[j] - vertices[i];

            // The normal of the edge depends on the polygon winding of the simplex
            norm = (vec2(line.y, -line.x) * (winding == CLOCKWISE ? 1 : -1)).unit();

            // Only keep the edge closest to the origin
            f64 dist = norm.dot(vertices[i]);
            if (dist < minDistance) {
                minDistance = dist;
                minNormal = norm;
                minIndex = j;
            }
        }
        return SimplexEdge(minDistance, minIndex, minNormal);
    }

    /** Returns true if the shapes collide */
    std::optional<GJKResult> getIntersection(Collider *colliderA, Collider *colliderB) {

        f64 e0 = (vertices[1].x - vertices[0].x)
            * (vertices[1].y + vertices[0].y);
        f64 e1 = (vertices[2].x - vertices[1].x)
            * (vertices[2].y + vertices[1].y);
        f64 e2 = (vertices[0].x - vertices[2].x)
            * (vertices[0].y + vertices[2].y);
        auto winding = (e0 + e1 + e2 >= 0)
            ? CLOCKWISE
            : COUNTERCLOCKWISE;

        vec2 intersection;
        for (int i = 0; i < MAX_TRIES; i++) {
            SimplexEdge edge = findClosestEdge(winding);
            // Calculate the difference for the two vertices furthest along the direction of the edge normal
            vec2 support = calculateSupport(colliderA, colliderB, edge.normal);
            // Check distance to the origin
            f64 distance = support.dot(edge.normal);
            intersection = edge.normal * distance;

            // If close enough, return if we need to move a distance greater than 0
            if (abs(distance - edge.distance) <= EPSILON) {
                f64 len = intersection.length();
                if (len != 0) return GJKResult(-intersection/len, len);
                else return {};
            }
            else {
                if (edge.index == 2) vertices[2] = support;
                else {
                    vertices[2] = vertices[1];
                    if (edge.index == 1) vertices[1] = support;
                    else {
                        vertices[1] = vertices[0];
                        vertices[0] = support;
                    }
                }
                length++;
            }
        }

        // Return if we need to move a distance greater than 0
        // Since we did more than the maximum amount of iterations, this may not be optimal
        f64 len = intersection.length();
        if (len != 0) return GJKResult(-intersection/len, len);
        else return {};
    }

    /** Returns a collision result if there was a collision */
    std::optional<GJKResult> gjkIntersection(Collider *colliderA, Collider *colliderB ) {
        vec2 direction = colliderB->center - colliderA->center;
        auto result = EVOLVE_IN_PROGRESS;
        while (result == EVOLVE_IN_PROGRESS) result = evolve(colliderA, colliderB, &direction);
        if (result != EVOLVE_DONE_FOUND_INTERSECTION) return {};
        return getIntersection(colliderA, colliderB);
    }
};

// MARK: Sweep-and-prune

class GameObj;
class SAPEdge {
    public:
    GameObj *object;
    bool isLeft;
    f64 x;
    SAPEdge(GameObj *object, bool isLeft, f64 x): object(object), isLeft(isLeft), x(x) {}
};


class GameObj {
    public:
    static GameObj *recentGameObjPull;
    int id;
    Collider *localArea;
    Collider *worldArea;
    SAPEdge *left;
    SAPEdge *right;
    vec2 areaOffset;
    vec2 areaScale;
    mat23 transform;
    GameObj(int id, SAPEdge *left, SAPEdge *right, Collider *localArea):
        id(id), localArea(localArea), worldArea(NULL), left(left), right(right),
        areaOffset(vec2(0, 0)), areaScale(vec2(1, 1)), transform(mat23()) {}
    void update() {
        GameObj::recentGameObjPull = this;
        pullTransform(id);
        mat23 t2 = (transform + areaOffset) * areaScale;
        localArea->transform(t2, &worldArea);
        EdgeBounds b = worldArea->bounds();
        left->x = b.left;
        right->x = b.right;
    }
};

GameObj *GameObj::recentGameObjPull = NULL;

typedef struct {
    List *touching;
    Simplex *simp;
    GameObj *obj;
} iterateStuff;

f64 gameobj_same_id_comparator(void *id, void *obj) { return *(int *)id - ((GameObj *)obj)->id; }

class SweepAndPrune {
    public:
    List edges;
    List objects;
    SweepAndPrune() {};
    void add(int id, Collider *coll) {
        GameObj *existing = objects.find(&id, gameobj_same_id_comparator);
        if (existing != NULL) {
            delete existing->localArea;
            delete existing->worldArea;
            existing->localArea = coll;
            existing->worldArea = NULL;
            return;
        }
        SAPEdge *left = new SAPEdge(NULL, true, NAN);
        SAPEdge *right = new SAPEdge(NULL, false, NAN);
        GameObj *obj = new GameObj(id, left, right, coll);
        left->object = right->object = obj;
        edges.push((void *)left);
        edges.push((void *)right);
        objects.push((void *)obj);
    }
    void remove(int id) {
        GameObj *obj = (GameObj *)objects.remove(&id, gameobj_same_id_comparator);
        if (obj != NULL) {
            // must call twice to remove both
            comparator hasEdge = [](void *obj, void *edge) -> f64 { return 1 - ((SAPEdge *)edge == ((GameObj *)obj)->left || (SAPEdge *)edge == ((GameObj *)obj)->right); };
            edges.remove((void *)obj, hasEdge);
            edges.remove((void *)obj, hasEdge);
            delete obj->left;
            delete obj->right;
            delete obj;
        }
    }
    void clear() {
        edges.clear();
        objects.clear();
    }
    void run() {
        // update all of the collider and edge data
        objects.iterate(NULL, [](void *arg, void *object) -> void {
            (void)arg;
            if (isObjActive(((GameObj *)object)->id))
                ((GameObj *)object)->update();
        });
        // sort edges
        edges.insertionSort([](void *a, void *b) -> f64 { return ((SAPEdge *)b)->x - ((SAPEdge *)a)->x; });
        // then find all the pairs that may be colliding
        List touching;
        Simplex simplex;
        iterateStuff stuff = {&touching, &simplex, NULL};
        edges.iterate((void *)&stuff, [](void *arg, void *item) -> void {
            SAPEdge *edge = (SAPEdge *)item;
            iterateStuff *stuff = (iterateStuff *)arg;
            GameObj *object = edge->object;
            if (edge->isLeft) {
                if (isObjActive(object->id)) {
                    stuff->obj = object;
                    stuff->touching->iterate(stuff, [](void *arg, void *item) -> void {
                        GameObj *other = (GameObj *)item;
                        iterateStuff *stuff = (iterateStuff *)arg;
                        if (isObjActive(stuff->obj->id) && !checkCollisionIgnore(other->id, stuff->obj->id)) {
                            // we found a colliding pair
                            stuff->simp->reset();
                            std::optional<GJKResult> res = stuff->simp->gjkIntersection(other->worldArea, stuff->obj->worldArea);
                            if (res.has_value())
                                handleCollision(other->id, stuff->obj->id, res->normal.x, res->normal.y, res->distance);
                        }
                    });
                }
                stuff->touching->push((void *)object);
            } else {
                stuff->touching->remove((void *)object, [](void *a, void *b) -> f64 { return (int)a - (int)b; });
            }
        });
    }
};

// MARK: MAIN EXPORTS

// global SAP instance
SweepAndPrune sap;

EXPORT void remove(int id) {
    sap.remove(id);
}

EXPORT void runCollision() {
    sap.run();
}

EXPORT void clear() {
    sap.clear();
}

EXPORT Collider *collider_allocate(ColliderType type, f64 cx, f64 cy) {
    switch (type) {
        case ELLIPSE: return new Ellipse(vec2(cx, cy), vec2(NAN, NAN), NAN);
        case POLYGON: return new Polygon(vec2(cx, cy));
        default: abort();
    }
}

EXPORT void registerObject(int objID, Collider *coll) {
    sap.add(objID, coll);
}

EXPORT void ellipse_sendData(Ellipse *e, f64 a, f64 rx, f64 ry) {
    e->angle = a;
    e->radii = vec2(rx, ry);
}

EXPORT vec2 *polygon_ensureVertices(Polygon *p, int nverts, bool isRect) {
    p->ensureCapacity(nverts);
    p->isRect = isRect;
    return p->vertices;
}

EXPORT void sendTransform(
    int id,
    f64 ox, f64 oy,
    f64 sx, f64 sy,
    f64 ax, f64 ay,
    f64 ma, f64 mb, f64 mc, f64 md, f64 me, f64 mf) {
        GameObj *target = GameObj::recentGameObjPull != NULL && GameObj::recentGameObjPull->id == id
                            ? GameObj::recentGameObjPull
                            : (GameObj *)sap.objects.find(&id, gameobj_same_id_comparator);
        if (target == NULL) abort();
        if (target->localArea == NULL) abort();
        target->localArea->anchor = vec2(ax, ay);
        target->areaOffset = vec2(ox, oy);
        target->areaScale = vec2(sx, sy);
        target->transform = mat23(ma, mb, mc, md, me, mf);
}
