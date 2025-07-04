import * as esbuild from "esbuild";
import * as fs from "node:fs";
import * as path from "node:path/posix";
import * as subprocess from "child_process"

// check for emscripten
try {
    subprocess.execSync("emcc --version >/dev/null");
} catch (e) {
    throw new Error(
        `

Emscripten is not installed but is required to compile!
Please go to https://emscripten.org/docs/getting_started/downloads.html
and install Emscripten.

`);
}

/** @type {esbuild.BuildOptions} */
const config = {
    bundle: true,
    sourcemap: true,
    minify: true,
    keepNames: false,
    metafile: true,
    platform: "browser",
    charset: "utf8",
    entryPoints: ["./index.ts"],
    format: "esm",
    target: "esnext",
    treeShaking: true,
    outfile: "build/plugin.js",
    plugins: [
        {
            name: "emscripten_wasm",
            setup(build) {
                build.onLoad({ filter: /\.cpp$/ }, async args => {
                    const wasmPath = args.path.replace(/\.cpp$/, ".wasm");
                    try {
                        subprocess.execFileSync("emcc", [
                            "-msimd128", // enable SIMD v128 instructions
                            "-std=c++23", "-stdlib=libc++", // use c++
                            "-Oz", // ultra optimizations
                            ...(process.argv.includes("-d") ? ["-g"] : []), // use debug symbols or not
                            "-sMALLOC=emmalloc", // tiny wasm-optimized malloc
                            "-sINITIAL_MEMORY=131072", // minimum memory to start with (2 pages, not 1, because of static vars for some reason)
                            "-sALLOW_MEMORY_GROWTH=1", // we're using malloc
                            "-sSTANDALONE_WASM=1", // don't import memory
                            "--no-entry", // don't expect a main()
                            "-flto", // turn on link-time optimization (smaller etc.)
                            "-o", wasmPath,
                            args.path
                        ]);
                    } catch (e) {
                        console.log("\n\n");
                        throw new Error(`compiling ${path.basename(args.path)} failed (exit code ${e.status})`);
                    }
                    console.info(`compiled ${path.basename(wasmPath)}<--${path.basename(args.path)}`);
                    return {
                        contents: new Uint8Array(fs.readFileSync(wasmPath).buffer),
                        loader: "binary"
                    }
                });
            }
        },
    ]
};

if (process.argv.includes("-w")) {
    config.plugins.push({
        name: "logger",
        setup(build) {
            build.onEnd(result => {
                if (result.errors.length == 0)
                    console.info(`[${new Date().toLocaleTimeString()}] rebuilt ${config.outfile} OK`);
                else
                    console.error(`[${new Date().toLocaleTimeString()}] failed to build ${config.outfile} !`);
            });
        },
    });
    await esbuild.context(config).then(ctx => ctx.watch());
}
else {
    const result = await esbuild.build(config);
    if (result.metafile) fs.writeFileSync("build/meta.json", JSON.stringify(result.metafile))
}
