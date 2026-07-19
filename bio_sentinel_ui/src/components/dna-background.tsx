import { useMemo } from "react";

function seededUnit(seed: number): number {
  const value = Math.sin(seed * 12.9898 + 78.233) * 43758.5453;
  return value - Math.floor(value);
}

/**
 * Ambient DNA stream background — subtle vertical columns of A/T/C/G
 * drifting slowly upward. Purely decorative; aria-hidden.
 */
export function DnaBackground({ columns = 28 }: { columns?: number }) {
  const streams = useMemo(() => {
    const letters = ["A", "T", "C", "G"];
    return Array.from({ length: columns }).map((_, i) => {
      const chars = Array.from({ length: 60 }).map(
        (_, index) => letters[Math.floor(seededUnit(i * 101 + index + 1) * 4)],
      );
      // duplicate for seamless loop
      const seq = [...chars, ...chars];
      return {
        id: i,
        left: (i / columns) * 100 + seededUnit(i * 17 + 2) * (100 / columns) * 0.4,
        duration: 60 + seededUnit(i * 17 + 3) * 90,
        delay: -seededUnit(i * 17 + 4) * 60,
        opacity: 0.04 + seededUnit(i * 17 + 5) * 0.02,
        fontSize: 12 + seededUnit(i * 17 + 6) * 4,
        chars: seq,
      };
    });
  }, [columns]);

  return (
    <div
      aria-hidden
      className="pointer-events-none fixed inset-0 overflow-hidden select-none"
      style={{ zIndex: 0 }}
    >
      {streams.map((s) => (
        <div
          key={s.id}
          className="dna-stream absolute top-0 flex flex-col items-center leading-[1.4] font-mono-tabular"
          style={{
            left: `${s.left}%`,
            animationDuration: `${s.duration}s`,
            animationDelay: `${s.delay}s`,
            color: "var(--color-dna)",
            opacity: s.opacity,
            fontSize: `${s.fontSize}px`,
            letterSpacing: "0.15em",
          }}
        >
          {s.chars.map((c, idx) => (
            <span key={idx}>{c}</span>
          ))}
        </div>
      ))}
    </div>
  );
}
