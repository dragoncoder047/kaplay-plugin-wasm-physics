import type { AreaComp, AreaCompOpt, GameObj, KAPLAYCtx, Vec2 } from "kaplay";
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
    readonly CT_CIRCLE: WebAssembly.Global<"i32">;
    readonly CT_ELLIPSE: WebAssembly.Global<"i32">;
    readonly CT_POLYGON: WebAssembly.Global<"i32">;
    runCollision(): void;
    /** The object takes ownership of the collider pointer */
    registerObject(id: number, colliderPtr: number): void;
    remove(id: number): void;
    clear(): void;

    /** returns a new thing, initialized all other fields to NaN. call another to update it */
    collider_allocate(type: number, x: number, y: number): number;
    circle_sendData(objPtr: number, radius: number): void;
    ellipse_sendData(objPtr: number, angle: number, radiusX: number, radiusY: number): void;
    /** returns pointer to the vertices which is a Float64Array with x, y, x, y, x, y, x, y, etc.
     * the first element (center) was already initialized. */
    polygon_ensureVertices(objPtr: number, numVertices: number): number
}

function isPaused(obj: GameObj): boolean {
    if (obj.paused) return true;
    return obj.parent ? isPaused(obj.parent) : false;
}

export function kaplayPhysicsWasm(k: KAPLAYCtx): Partial<KAPLAYCtx> {
    var wasm: PhysicsWasmExports;
    const idToObjMap = new Map<number, GameObj<AreaComp>>();


    class Collision {
        source: GameObj;
        target: GameObj;
        normal: Vec2;
        distance: number;
        resolved: boolean = false;
        constructor(
            source: GameObj,
            target: GameObj,
            normal: Vec2,
            distance: number,
            resolved = false,
        ) {
            this.source = source;
            this.target = target;
            this.normal = normal;
            this.distance = distance;
            this.resolved = resolved;
        }
        get displacement() {
            return this.normal.scale(this.distance);
        }
        reverse() {
            return new Collision(
                this.target,
                this.source,
                this.normal.scale(-1),
                this.distance,
                this.resolved,
            );
        }
        hasOverlap() {
            return this.distance > 0;
        }
        isLeft() {
            return this.normal.cross(k._k.game.gravity || k.vec2(0, 1)) > 0;
        }
        isRight() {
            return this.normal.cross(k._k.game.gravity || k.vec2(0, 1)) < 0;
        }
        isTop() {
            return this.normal.dot(k._k.game.gravity || k.vec2(0, 1)) > 0;
        }
        isBottom() {
            return this.normal.dot(k._k.game.gravity || k.vec2(0, 1)) < 0;
        }
        preventResolution() {
            this.resolved = true;
        }
    }



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
                const col1 = new Collision(
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
            emscripten_notify_memory_growth() {
                k.debug.log("WASM module called (memory.grow) successfully")
            }
        }
    };



    k.load((async () => {
        wasm = await initWASM(wasmAPIFuncs);
        wasm._initialize();
    })());

    k.system("collision", () => {
        wasm.runCollision();
    }, [k.SystemPhase.AfterFixedUpdate]);

    k.onSceneLeave(() => {
        wasm.clear();
    });

    const oldArea = k.area;
    return {
        area(opt?: AreaCompOpt) {
            const comp = oldArea(opt);
            const { add: oldAdd, destroy: oldDestroy } = comp;
            Object.assign(comp, {
                add(this: GameObj<AreaComp>) {
                    oldAdd?.call(this);
                    idToObjMap.set(this.id, this);
                    wasm.registerObject(this.id);
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
