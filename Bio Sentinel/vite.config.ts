// @lovable.dev/vite-tanstack-config already includes TanStack Start, React,
// Tailwind, path aliases, and Nitro. Keep the integration config here narrow:
// only proxy API traffic during local development.
import { defineConfig } from "@lovable.dev/vite-tanstack-config";
import process from "node:process";

const genomeFirewallApiUrl = process.env.GENOME_FIREWALL_API_URL || "http://127.0.0.1:8000";

export default defineConfig({
  vite: {
    server: {
      proxy: {
        "/api": {
          target: genomeFirewallApiUrl,
          changeOrigin: true,
        },
      },
    },
  },
  tanstackStart: {
    server: { entry: "server" },
  },
});
