#include <emscripten.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <optional>

#define EXPORT EMSCRIPTEN_KEEPALIVE extern "C"

// MARK: imports

/** handles the collisions, when found */
EM_IMPORT(handleCollision) void handleCollision(int idA, int idB, double normX, double normY, double distance);

/** returns true if the objects shouldn't collide with each other. */
EM_IMPORT(checkCollisionIgnore) bool checkCollisionIgnore(int idA, int idB);

/** returns true if the object isn't paused */
EM_IMPORT(isObjActive) bool isObjActive(int id);

/* tells the JS code to send the transformation data for the object with this id */
EM_IMPORT(pullTransform) void pullTransform(int id);

// MARK: lists n' stuff

typedef void (*iter_body)(void *arg, void *item);
typedef double (*comparator)(void *left, void *right);

class List {
    public:
    class Node {
        public:
        void *item;
        Node *next;
        Node(void *item, Node *next): item(item), next(next) {}
    };
    Node *head;
    List(): head(NULL) {}
    ~List() { clear(); }
    void clear() {
        while (head != NULL) {
            Node *next = head->next;
            delete head;
            head = next;
        }
    }
    void push(void *item) { head = new Node(item, head); }
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
        Node *current = head;
        while (current != NULL) {
            Node *next = current->next;
            if (sorted == NULL || compare(sorted->item, current->item) < 0) {
                current->next = sorted;
                sorted = current;
            } else {
                Node *after = sorted;
                while (after->next != NULL && compare(after->next->item, current->item) >= 0)
                    after = after->next;
                current->next = after->next;
                after->next = current;
            }
            current = next;
        }
        head = sorted;
    }
};

// MARK: mathy math math

class vec2 {
    public:
    double x, y;
    vec2(): x(NAN), y(NAN) {}
    vec2(double x, double y): x(x), y(y) {}
    double length() const { return sqrt(x*x + y*y); }
    inline vec2 operator+(vec2 other) const { return vec2(x + other.x, y + other.y); }
    inline vec2 operator-(vec2 other) const { return vec2(x - other.x, y - other.y); }
    inline vec2 operator-() const { return vec2(-x, -y); }
    inline vec2 operator*(double other) const { return vec2(x * other, y * other); }
    inline vec2 operator*(vec2 other) const { return vec2(x * other.x, y * other.y); }
    inline vec2 operator/(double other) const { return vec2(x / other, y / other); }
    inline vec2 unit() const { double l = length(); return vec2(x/l, y/l); }
    inline double dot(vec2 other) const { return x * other.x + y * other.y; }
    inline vec2 rotate(double angle) const { double s = sin(angle), c = cos(angle); return vec2(x*c - y*s, x*s + y*c); }
};

class mat23 {
    public:
    double a, b, c, d, e, f;
    mat23(): a(1), b(0), c(0), d(1), e(0), f(0) {}
    mat23(double a, double b, double c, double d, double e, double f):
        a(a), b(b), c(c), d(d), e(e), f(f) {}
    inline mat23 operator+(vec2 t) const { return mat23(a, b, c, d, e + vec2(a, c).dot(t), f + vec2(b, d).dot(t)); }
    inline mat23 operator*(vec2 s) const { return mat23(a * s.x, b * s.x, c * s.y, d * s.y, e, f); }
    inline vec2 transformVector(vec2 v) const { return vec2(vec2(a, c).dot(v), vec2(b, d).dot(v)); }
    inline vec2 transformPoint(vec2 v) const { return vec2(e, f) + transformVector(v); }
};

// MARK: GJK algorithm

struct XBounds { double left, right; };

enum ColliderType {
    NONE,
    CIRCLE,
    ELLIPSE,
    POLYGON
};
EXPORT ColliderType CT_NONE = ColliderType::NONE;
EXPORT ColliderType CT_CIRCLE = ColliderType::CIRCLE;
EXPORT ColliderType CT_ELLIPSE = ColliderType::ELLIPSE;
EXPORT ColliderType CT_POLYGON = ColliderType::POLYGON;

class Collider {
    public:
    vec2 center;
    ColliderType type;
    Collider(): type(NONE) {}
    Collider(vec2 center, ColliderType type): center(center), type(type) {}
    virtual vec2 support(vec2 direction) const { return vec2(NAN, NAN); };
    virtual XBounds bounds() const { return {NAN, NAN}; }
    virtual void transform(Collider **out) { abort(); }
};

class Circle: public Collider {
    public:
    double radius;
    Circle(vec2 center, double radius): Collider(center, CIRCLE), radius(radius) { }
    vec2 support(vec2 direction) const { return direction.unit() * radius + center; }
    XBounds bounds() const { return {center.x - radius, center.x + radius}; }
};

class Ellipse: public Collider {
    public:
    double angle;
    vec2 radii;
    Ellipse(vec2 center, vec2 radii, double angle):
        Collider(center, ELLIPSE), radii(radii), angle(angle) {}
    vec2 support(vec2 direction) {
        return (center
            + (angle == 0
                // Axis aligned
                ? (direction.unit() * radii)
                // Rotated
                : (direction.rotate(-angle).unit() * radii).rotate(angle)));
    }
    XBounds bounds() {
        if (angle == 0) {
            return {center.x - radii.x, center.x + radii.x};
        }
        double halfwidth = (radii * vec2(1, 0).rotate(angle)).length();
        return {center.x - halfwidth, center.x + halfwidth};
    }
};

class Polygon: public Collider {
    public:
    vec2 *vertices;
    int length = 0, capacity = 0;
    Polygon(vec2 center): Collider(center, POLYGON) { addVertex(center); }
    ~Polygon() { free(vertices); }
    void ensureCapacity(int length) {
        if (length <= capacity) return;
        while (length > capacity) capacity <<= 1;
        vec2 *newVertices = (vec2 *)calloc(capacity, sizeof(vec2));
        memcpy(newVertices, vertices, length * sizeof(vec2));
        free(vertices);
        vertices = newVertices;
    }
    void addVertex(vec2 v) {
        ensureCapacity(length + 1);
        vertices[length++] = v;
    }
    vec2 support(vec2 direction) {
        vec2 maxPoint = vec2(NAN, NAN);
        double maxDistance = -INFINITY;
        for (size_t i = 0; i < length; i++) {
            vec2 vertex = vertices[i];
            double dist = vertex.dot(direction);
            if (dist > maxDistance) {
                maxDistance = dist;
                maxPoint = vertex;
            }
        }
        return maxPoint;
    }
    XBounds bounds() {
        double min = INFINITY, max = -INFINITY;
        for (size_t i = 0; i < length; i++) {
            vec2 vertex = vertices[i];
            if (vertex.x > max) max = vertex.x;
            if (vertex.x < min) min = vertex.x;
        }
        return {min, max};
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
    double n = a.x * b.y - a.y * b.x;
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
    double distance;
    int index;
    vec2 normal;
    SimplexEdge(double distance, int index, vec2 normal): distance(distance), index(index), normal(normal) {}
};

class GJKResult {
    public:
    vec2 normal;
    double distance;
    GJKResult(vec2 normal, double distance): normal(normal), distance(distance) {}
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
        double minDistance = INFINITY;
        vec2 minNormal, line, norm;
        int minIndex = 0;
        for (int i = 0; i < length; i++) {
            int j = (i + 1) % length;
            line = vertices[j] - vertices[i];

            // The normal of the edge depends on the polygon winding of the simplex
            norm = (vec2(line.y, -line.x) * (winding == CLOCKWISE ? 1 : -1)).unit();

            // Only keep the edge closest to the origin
            double dist = norm.dot(vertices[i]);
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
        const double EPSILON = 0.00001;

        double e0 = (vertices[1].x - vertices[0].x)
            * (vertices[1].y + vertices[0].y);
        double e1 = (vertices[2].x - vertices[1].x)
            * (vertices[2].y + vertices[1].y);
        double e2 = (vertices[0].x - vertices[2].x)
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
            double distance = support.dot(edge.normal);
            intersection = edge.normal * distance;

            // If close enough, return if we need to move a distance greater than 0
            if (abs(distance - edge.distance) <= EPSILON) {
                double len = intersection.length();
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
        double len = intersection.length();
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
    double x;
    SAPEdge(GameObj *object, bool isLeft, double x): object(object), isLeft(isLeft), x(x) {}
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
        id(id), left(left), right(right), localArea(localArea), worldArea(NULL),
        areaOffset(vec2(0, 0)), areaScale(vec2(1, 1)), transform(mat23()) {}
    void updateEdges() {
        if (!isObjActive(id)) return;
        GameObj::recentGameObjPull = this;
        pullTransform(id);
        // TODO: update world area transform
        abort();
    }
};

GameObj *GameObj::recentGameObjPull = NULL;

typedef struct {
    List *touching;
    Simplex *simp;
    GameObj *obj;
} iterateStuff;

double gameobj_same_id_comparator(void *id, void *obj) { return *(int *)id - ((GameObj *)obj)->id; }

class SweepAndPrune {
    public:
    List edges;
    List objects;
    SweepAndPrune() {};
    void add(int id, Collider *coll) {
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
            comparator hasEdge = [](void *obj, void *edge) -> double { return 1 - ((SAPEdge *)edge == ((GameObj *)obj)->left || (SAPEdge *)edge == ((GameObj *)obj)->right); };
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
        objects.iterate(NULL, [](void *arg, void *object) -> void { (void)arg; ((GameObj *)object)->updateEdges(); });
        // sort edges
        edges.insertionSort([](void *a, void *b) -> double { return ((SAPEdge *)b)->x - ((SAPEdge *)a)->x; });
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
                stuff->touching->remove((void *)object, [](void *a, void *b) -> double { return (int)a - (int)b; });
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

EXPORT Collider *collider_allocate(ColliderType type, double cx, double cy) {
    switch (type) {
        case CIRCLE: return new Circle(vec2(cx, cy), NAN);
        case ELLIPSE: return new Ellipse(vec2(cx, cy), vec2(NAN, NAN), NAN);
        case POLYGON: return new Polygon(vec2(cx, cy));
        default: abort();
    }
}

EXPORT void registerObject(int objID, Collider *coll) {
    sap.add(objID, coll);
}

EXPORT void circle_sendData(Circle *c, double r) {
    c->radius = r;
}

EXPORT void ellipse_sendData(Ellipse *e, double a, double rx, double ry) {
    e->angle = a;
    e->radii = vec2(rx, ry);
}

EXPORT vec2 *polygon_ensureVertices(Polygon *p, int nverts) {
    p->ensureCapacity(nverts);
    return p->vertices;
}

EXPORT void sendTransform(
    int id,
    double ox, double oy,
    double sx, double sy,
    double ma, double mb, double mc, double md, double me, double mf) {
        GameObj *target = GameObj::recentGameObjPull != NULL && GameObj::recentGameObjPull->id == id
                            ? GameObj::recentGameObjPull
                            : (GameObj *)sap.objects.find(&id, gameobj_same_id_comparator);
        if (target == NULL) abort();
        target->areaOffset = vec2(ox, oy);
        target->areaScale = vec2(sx, sy);
        target->transform = mat23(ma, mb, mc, md, me, mf);
}
