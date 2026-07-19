import { useMemo } from "react";

/**
 * Ambient DNA stream background — subtle vertical columns of A/T/C/G
 * drifting slowly upward. Purely decorative; aria-hidden.
 */
export function DnaBackground({ columns = 28 }: { columns?: number }) {
  const streams = useMemo(() => {
    const letters = ["A", "T", "C", "G"];
    return Array.from({ length: columns }).map((_, i) => {
      const chars = Array.from({ length: 60 }).map(
        () => letters[Math.floor(Math.random() * 4)],
      );
      // duplicate for seamless loop
      const seq = [...chars, ...chars];
      return {
        id: i,
        left: (i / columns) * 100 + Math.random() * (100 / columns) * 0.4,
        duration: 60 + Math.random() * 90,
        delay: -Math.random() * 60,
        opacity: 0.04 + Math.random() * 0.02,
        fontSize: 12 + Math.random() * 4,
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
