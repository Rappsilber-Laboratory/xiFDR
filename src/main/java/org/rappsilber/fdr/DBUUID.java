package org.rappsilber.fdr;
import java.util.UUID;
import java.util.concurrent.atomic.AtomicLong;
import java.security.SecureRandom;

public class DBUUID {
    private static final long PREFIX = 0xFFFFFFFF00000000L;
    private static final SecureRandom RANDOM = new SecureRandom();
    private static final AtomicLong COUNTER = new AtomicLong(0);

    public static UUID dbUUID() {
        // Upper 64 bits: prefix + timestamp
        long timestamp = System.currentTimeMillis() / 1000;  // seconds since epoch
        long uuidUpper = PREFIX | (timestamp & 0xFFFFFFFFL);

        // Lower 64 bits: counter (top 8 bits) + random (lower 56 bits)
        long counter = COUNTER.getAndIncrement() % 64;
        long randomValue = RANDOM.nextLong() & 0x00FFFFFFFFFFFFFFL; // 56 random bits
        long uuidLower = (counter << 56) | randomValue;

        // Combine into a UUID
        return new UUID(uuidUpper, uuidLower);
    }
}