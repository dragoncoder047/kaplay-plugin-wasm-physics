import type { AnchorComp, AreaComp, AreaCompOpt, GameObj, KAPLAYCtx, Shape, Vec2 } from "kaplay";
import physicsWasmBin from "./physics.cpp";

async function initWASM(imports: Record<string, any>): Promise<PhysicsWasmExports> {
    const mod = await WebAssembly.compile(physicsWasmBin);
    const inst = await WebAssembly.instantiate(mod, imports);
    return inst.exports as PhysicsWasmExports;
}

interface PhysicsWasmExports extends WebAssembly.Exports {
    /** runs all of the C++ top-level constructors */
    _initialize(): void;
    readonly memory: WebAssembly.Memory;
    readonly CT_ELLIPSE: WebAssembly.Global<"i32">;
    readonly CT_POLYGON: WebAssembly.Global<"i32">;
    runCollision(): void;
    /** The object takes ownership of the collider pointer */
    registerObject(id: number, colliderPtr: number): void;
    remove(id: number): void;
    clear(): void;

    /** returns a new thing, initialized all other fields to NaN. call another to update it */
    collider_allocate(type: number, x: number, y: number): number;
    ellipse_sendData(objPtr: number, angle: number, radiusX: number, radiusY: number): void;
    /** returns pointer to the vertices which is a Float64Array with x, y, x, y, x, y, x, y, etc.
     * the first element (center) was already initialized. */
    polygon_ensureVertices(objPtr: number, numVertices: number, isRect: number): number

    sendTransform(
        id: number,
        ox: number, oy: number,
        sx: number, sy: number,
        ax: number, ay: number,
        ma: number, mb: number, mc: number, md: number, me: number, mf: number): void;
}

function isPaused(obj: GameObj): boolean {
    if (obj.paused) return true;
    return obj.parent ? isPaused(obj.parent) : false;
}

export function kaplayPhysicsWasm(k: KAPLAYCtx): Partial<KAPLAYCtx> {
    var wasm: PhysicsWasmExports;
    const idToObjMap = new Map<number, GameObj<AreaComp | AnchorComp>>();

    const wasmAPIFuncs: WebAssembly.Imports = {
        env: { // must be called env because EM_IMPORT hardcodes the module name to "env"
            isObjActive(id: number) {
                const obj = idToObjMap.get(id);
                if (!obj) return 0; // wtf?? skip it
                return obj.exists() && !isPaused(obj) ? 1 : 0;
            },
            checkCollisionIgnore(idA: number, idB: number) {
                const obj = idToObjMap.get(idA)!;
                const other = idToObjMap.get(idB)!;
                for (const tag of obj.collisionIgnore) if (other.is(tag)) return 1;
                for (const tag of other.collisionIgnore) if (obj.is(tag)) return 1;
                return 0;
            },
            handleCollision(idA: number, idB: number, normX: number, normY: number, distance: number) {
                const obj = idToObjMap.get(idA)!;
                const other = idToObjMap.get(idB)!;
                const col1 = new k.Collision(
                    obj,
                    other,
                    k.vec2(normX, normY),
                    distance,
                );
                obj.trigger("collideUpdate", other, col1);
                const col2 = col1.reverse();
                // resolution only has to happen once
                col2.resolved = col1.resolved;
                other.trigger("collideUpdate", obj, col2);
            },
            pullTransform(id: number) {
                const o = idToObjMap.get(id);
                if (!o) throw new Error("object not found!!");
                const { x: sx, y: sy } = o.area.scale;
                const { x: ox, y: oy } = o.area.offset;
                const { x: ax, y: ay } = k.anchorToVec2(o.anchor ?? "topleft");
                const { a, b, c, d, e, f } = o.transform;
                wasm.sendTransform(id, ox, oy, sx, sy, ax, ay, a, b, c, d, e, f);
            },
            emscripten_notify_memory_growth(index: number) {
                k.debug.log(`WASM module called (memory.grow) for memory #${index}`)
            }
        }
    };

    /**
     * WARNING: MUST pass the return value of this directly to wasm.registerObject() or there will be a memory leak!!!!!!!
     */
    const allocWasmAreaObj = (localArea: Shape): number => {
        if (localArea instanceof k.Circle) {
            const circlePtr = wasm.collider_allocate(wasm.CT_ELLIPSE.value, localArea.center.x, localArea.center.y);
            wasm.ellipse_sendData(circlePtr, 0, localArea.radius, localArea.radius);
            return circlePtr;
        } else if (localArea instanceof k.Ellipse) {
            const ellipsePtr = wasm.collider_allocate(wasm.CT_ELLIPSE.value, localArea.center.x, localArea.center.y);
            wasm.ellipse_sendData(ellipsePtr, k.deg2rad(localArea.angle), localArea.radiusX, localArea.radiusY);
            return ellipsePtr;
        } else {
            const points: [Vec2, ...Vec2[]] = (
                localArea instanceof k.Rect ? localArea.points() :
                    localArea instanceof k.Polygon ? localArea.pts as any :
                        localArea instanceof k.Line ? [localArea.p1, localArea.p1, localArea.p2] :
                            localArea instanceof k.Point ? [localArea.pt, localArea.pt, localArea.pt] :
                                (() => { throw new Error("unknown shape type") })())
            const polygonPtr = wasm.collider_allocate(wasm.CT_POLYGON.value, points[0].x, points[0].y);
            const pointsArrPtr = wasm.polygon_ensureVertices(polygonPtr, points.length, +(localArea instanceof k.Rect));
            new Float64Array(wasm.memory.buffer, pointsArrPtr, points.length * 2)
                .set(points.flatMap(p => [p.x, p.y]));
            return polygonPtr;
        }
    };


    k.load((async () => {
        wasm = await initWASM(wasmAPIFuncs);
        wasm._initialize();
    })());
    k.system("collision", () => {
        wasm.runCollision();
    }, [k.SystemPhase.AfterFixedUpdate, k.SystemPhase.AfterUpdate]);

    const oldArea = k.area;
    return {
        area(opt?: AreaCompOpt) {
            const comp = oldArea(opt);
            const { add: oldAdd, destroy: oldDestroy } = comp;
            Object.assign(comp, {
                add(this: GameObj<AreaComp | AnchorComp>) {
                    oldAdd?.call(this);
                    idToObjMap.set(this.id, this);
                    wasm.registerObject(this.id, allocWasmAreaObj(this.localArea()));
                },
                destroy(this: GameObj<AreaComp>) {
                    oldDestroy?.call(this);
                    idToObjMap.delete(this.id);
                    wasm.remove(this.id);
                }
            });
            return comp;
        },
    }
}
export default kaplayPhysicsWasm;
