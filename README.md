# kaplay-plugin-wasm-physics

This is a [KAPLAY][] plugin that replaces the default area() system with one written in WebAssembly, which is much faster.

> [!WARNING]
> THIS IS CURRENTLY INCOMPLETE AND UNRELEASED! DO NOT INSTALL UNLESS YOU ARE DEVELOPING IT!!!!

## Installation

Not on NPM, so just install from GitHub directly using npm or [pnpm].

First make sure you have [Emscripten] installed so you can compile it.

```sh
pnpm add https://github.com/dragoncoder047/kaplay-plugin-wasm-physics
```

## Usage

Just install to KAPLAY like you would any other plugin:

```ts
import kaplayWasmPhysics from "kaplay-plugin-wasm-physics";

kaplay({
    // ...other options
    plugins: [
        kaplayWasmPhysics,
        // ...other plugins
    ],
});
```

It's a drop-in replacement, so no other API changes will occur.

[KAPLAY]: https://kaplayjs.com/
[Emscripten]: https://emscripten.org/
[pnpm]: https://pnpm.io/
