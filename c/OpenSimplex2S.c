/**
 * K.jpg's OpenSimplex 2, smooth variant ("SuperSimplex")
 */

#include "_OpenSimplex2S.h"

typedef int boolean;

#define FALSE 0
#define TRUE 1

static const float GRADIENTS_2D[] = {
    #include S_GRADIENTS_2D_SRC
};
static const float GRADIENTS_3D[] = {
    #include S_GRADIENTS_3D_SRC
};
static const float GRADIENTS_4D[] = {
    #include S_GRADIENTS_4D_SRC
};
static const int LOOKUP_4D_A[] = {
    #include S_LOOKUP_4D_A_SRC
};
static const LatticeVertex4D LOOKUP_4D_B[] = {
    #include S_LOOKUP_4D_B_SRC
};


/*
 * Utility
 */

static float noise2_grad(long seed, long xsvp, long ysvp, float dx, float dy) {
    long hash = seed ^ xsvp ^ ysvp;
    hash *= S_HASH_MULTIPLIER;
    hash ^= hash >> (64 - S_N_GRADS_2D_EXPONENT + 1);
    int gi = (int)hash & ((S_N_GRADS_2D - 1) << 1);
    return GRADIENTS_2D[gi | 0] * dx + GRADIENTS_2D[gi | 1] * dy;
}

static float noise3_grad(long seed, long xrvp, long yrvp, long zrvp, float dx, float dy, float dz) {
    long hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
    hash *= S_HASH_MULTIPLIER;
    hash ^= hash >> (64 - S_N_GRADS_3D_EXPONENT + 2);
    int gi = (int)hash & ((S_N_GRADS_3D - 1) << 2);
    return GRADIENTS_3D[gi | 0] * dx + GRADIENTS_3D[gi | 1] * dy + GRADIENTS_3D[gi | 2] * dz;
}

static float noise4_grad(long seed, long xsvp, long ysvp, long zsvp, long wsvp, float dx, float dy, float dz, float dw) {
    long hash = seed ^ (xsvp ^ ysvp) ^ (zsvp ^ wsvp);
    hash *= S_HASH_MULTIPLIER;
    hash ^= hash >> (64 - S_N_GRADS_4D_EXPONENT + 2);
    int gi = (int)hash & ((S_N_GRADS_4D - 1) << 2);
    return (GRADIENTS_4D[gi | 0] * dx + GRADIENTS_4D[gi | 1] * dy) + (GRADIENTS_4D[gi | 2] * dz + GRADIENTS_4D[gi | 3] * dw);
}

static int fastFloor(double x) {
    int xi = (int)x;
    return x < xi ? xi - 1 : xi;
}

/*
 * Noise Evaluators
 */

/**
 * 2D  OpenSimplex2S/SuperSimplex noise base.
 */
static float noise2_UnskewedBase(long seed, double xs, double ys) {

    // Get base points and offsets.
    int xsb = fastFloor(xs), ysb = fastFloor(ys);
    float xi = (float)(xs - xsb), yi = (float)(ys - ysb);

    // Prime pre-multiplication for hash.
    long xsbp = xsb * S_PRIME_X, ysbp = ysb * S_PRIME_Y;

    // Unskew.
    float t = (xi + yi) * (float)S_UNSKEW_2D;
    float dx0 = xi + t, dy0 = yi + t;

    // First vertex.
    float a0 = S_RSQUARED_2D - dx0 * dx0 - dy0 * dy0;
    float value = (a0 * a0) * (a0 * a0) * noise2_grad(seed, xsbp, ysbp, dx0, dy0);

    // Second vertex.
    float a1 = (float)(2 * (1 + 2 * S_UNSKEW_2D) * (1 / S_UNSKEW_2D + 2)) * t + ((float)(-2 * (1 + 2 * S_UNSKEW_2D) * (1 + 2 * S_UNSKEW_2D)) + a0);
    float dx1 = dx0 - (float)(1 + 2 * S_UNSKEW_2D);
    float dy1 = dy0 - (float)(1 + 2 * S_UNSKEW_2D);
    value += (a1 * a1) * (a1 * a1) * noise2_grad(seed, xsbp + S_PRIME_X, ysbp + S_PRIME_Y, dx1, dy1);

    // Third and fourth vertices.
    // Nested conditionals were faster than compact bit logic/arithmetic.
    float xmyi = xi - yi;
    if (t < S_UNSKEW_2D) {
        if (xi + xmyi > 1) {
            float dx2 = dx0 - (float)(3 * S_UNSKEW_2D + 2);
            float dy2 = dy0 - (float)(3 * S_UNSKEW_2D + 1);
            float a2 = S_RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0) {
                value += (a2 * a2) * (a2 * a2) * noise2_grad(seed, xsbp + (S_PRIME_X << 1), ysbp + S_PRIME_Y, dx2, dy2);
            }
        }
        else
        {
            float dx2 = dx0 - (float)S_UNSKEW_2D;
            float dy2 = dy0 - (float)(S_UNSKEW_2D + 1);
            float a2 = S_RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0) {
                value += (a2 * a2) * (a2 * a2) * noise2_grad(seed, xsbp, ysbp + S_PRIME_Y, dx2, dy2);
            }
        }

        if (yi - xmyi > 1) {
            float dx3 = dx0 - (float)(3 * S_UNSKEW_2D + 1);
            float dy3 = dy0 - (float)(3 * S_UNSKEW_2D + 2);
            float a3 = S_RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
            if (a3 > 0) {
                value += (a3 * a3) * (a3 * a3) * noise2_grad(seed, xsbp + S_PRIME_X, ysbp + (S_PRIME_Y << 1), dx3, dy3);
            }
        }
        else
        {
            float dx3 = dx0 - (float)(S_UNSKEW_2D + 1);
            float dy3 = dy0 - (float)S_UNSKEW_2D;
            float a3 = S_RSQUARED_2D - dx3 * dx3 - dy3 * dy3;
            if (a3 > 0) {
                value += (a3 * a3) * (a3 * a3) * noise2_grad(seed, xsbp + S_PRIME_X, ysbp, dx3, dy3);
            }
        }
    }
    else
    {
        if (xi + xmyi < 0) {
            float dx2 = dx0 + (float)(1 + S_UNSKEW_2D);
            float dy2 = dy0 + (float)S_UNSKEW_2D;
            float a2 = S_RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0) {
                value += (a2 * a2) * (a2 * a2) * noise2_grad(seed, xsbp - S_PRIME_X, ysbp, dx2, dy2);
            }
        }
        else
        {
            float dx2 = dx0 - (float)(S_UNSKEW_2D + 1);
            float dy2 = dy0 - (float)S_UNSKEW_2D;
            float a2 = S_RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0) {
                value += (a2 * a2) * (a2 * a2) * noise2_grad(seed, xsbp + S_PRIME_X, ysbp, dx2, dy2);
            }
        }

        if (yi < xmyi) {
            float dx2 = dx0 + (float)S_UNSKEW_2D;
            float dy2 = dy0 + (float)(S_UNSKEW_2D + 1);
            float a2 = S_RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0) {
                value += (a2 * a2) * (a2 * a2) * noise2_grad(seed, xsbp, ysbp - S_PRIME_Y, dx2, dy2);
            }
        }
        else
        {
            float dx2 = dx0 - (float)S_UNSKEW_2D;
            float dy2 = dy0 - (float)(S_UNSKEW_2D + 1);
            float a2 = S_RSQUARED_2D - dx2 * dx2 - dy2 * dy2;
            if (a2 > 0) {
                value += (a2 * a2) * (a2 * a2) * noise2_grad(seed, xsbp, ysbp + S_PRIME_Y, dx2, dy2);
            }
        }
    }

    return value;
}

/**
 * 2D OpenSimplex2S/SuperSimplex noise, standard lattice orientation.
 */
float OpenSimplex2S_noise2(long seed, double x, double y) {

    // Get points for A2* lattice
    double s = S_SKEW_2D * (x + y);
    double xs = x + s, ys = y + s;

    return noise2_UnskewedBase(seed, xs, ys);
}

/**
 * 2D OpenSimplex2S/SuperSimplex noise, with Y pointing down the main diagonal.
 * Might be better for a 2D sandbox style game, where Y is vertical.
 * Probably slightly less optimal for heightmaps or continent maps,
 * unless your map is centered around an equator. It's a slight
 * difference, but the option is here to make it easy.
 */
float OpenSimplex2S_noise2_ImproveX(long seed, double x, double y) {

    // Skew transform and rotation baked into one.
    double xx = x * S_ROOT2OVER2;
    double yy = y * (S_ROOT2OVER2 * (1 + 2 * S_SKEW_2D));

    return noise2_UnskewedBase(seed, yy + xx, yy - xx);
}

/**
 * Generate overlapping cubic lattices for 3D Re-oriented BCC noise.
 * Lookup table implementation inspired by DigitalShadow.
 * It was actually faster to narrow down the points in the loop itself,
 * than to build up the index with enough info to isolate 8 points.
 */
static float noise3_UnrotatedBase(long seed, double xr, double yr, double zr) {

    // Get base points and offsets.
    int xrb = fastFloor(xr), yrb = fastFloor(yr), zrb = fastFloor(zr);
    float xi = (float)(xr - xrb), yi = (float)(yr - yrb), zi = (float)(zr - zrb);

    // Prime pre-multiplication for hash. Also flip seed for second lattice copy.
    long xrbp = xrb * S_PRIME_X, yrbp = yrb * S_PRIME_Y, zrbp = zrb * S_PRIME_Z;
    long seed2 = seed ^ -0x52D547B2E96ED629L;

    // -1 if positive, 0 if negative.
    int xNMask = (int)(-0.5f - xi), yNMask = (int)(-0.5f - yi), zNMask = (int)(-0.5f - zi);

    // First vertex.
    float x0 = xi + xNMask;
    float y0 = yi + yNMask;
    float z0 = zi + zNMask;
    float a0 = S_RSQUARED_3D - x0 * x0 - y0 * y0 - z0 * z0;
    float value = (a0 * a0) * (a0 * a0) * noise3_grad(seed,
            xrbp + (xNMask & S_PRIME_X), yrbp + (yNMask & S_PRIME_Y), zrbp + (zNMask & S_PRIME_Z), x0, y0, z0);

    // Second vertex.
    float x1 = xi - 0.5f;
    float y1 = yi - 0.5f;
    float z1 = zi - 0.5f;
    float a1 = S_RSQUARED_3D - x1 * x1 - y1 * y1 - z1 * z1;
    value += (a1 * a1) * (a1 * a1) * noise3_grad(seed2,
            xrbp + S_PRIME_X, yrbp + S_PRIME_Y, zrbp + S_PRIME_Z, x1, y1, z1);

    // Shortcuts for building the remaining falloffs.
    // Derived by subtracting the polynomials with the offsets plugged in.
    float xAFlipMask0 = ((xNMask | 1) << 1) * x1;
    float yAFlipMask0 = ((yNMask | 1) << 1) * y1;
    float zAFlipMask0 = ((zNMask | 1) << 1) * z1;
    float xAFlipMask1 = (-2 - (xNMask << 2)) * x1 - 1.0f;
    float yAFlipMask1 = (-2 - (yNMask << 2)) * y1 - 1.0f;
    float zAFlipMask1 = (-2 - (zNMask << 2)) * z1 - 1.0f;

    boolean skip5 = FALSE;
    float a2 = xAFlipMask0 + a0;
    if (a2 > 0) {
        float x2 = x0 - (xNMask | 1);
        float y2 = y0;
        float z2 = z0;
        value += (a2 * a2) * (a2 * a2) * noise3_grad(seed,
                xrbp + (~xNMask & S_PRIME_X), yrbp + (yNMask & S_PRIME_Y), zrbp + (zNMask & S_PRIME_Z), x2, y2, z2);
    }
    else
    {
        float a3 = yAFlipMask0 + zAFlipMask0 + a0;
        if (a3 > 0) {
            float x3 = x0;
            float y3 = y0 - (yNMask | 1);
            float z3 = z0 - (zNMask | 1);
            value += (a3 * a3) * (a3 * a3) * noise3_grad(seed,
                    xrbp + (xNMask & S_PRIME_X), yrbp + (~yNMask & S_PRIME_Y), zrbp + (~zNMask & S_PRIME_Z), x3, y3, z3);
        }

        float a4 = xAFlipMask1 + a1;
        if (a4 > 0) {
            float x4 = (xNMask | 1) + x1;
            float y4 = y1;
            float z4 = z1;
            value += (a4 * a4) * (a4 * a4) * noise3_grad(seed2,
                    xrbp + (xNMask & ((unsigned long)S_PRIME_X * 2)), yrbp + S_PRIME_Y, zrbp + S_PRIME_Z, x4, y4, z4);
            skip5 = TRUE;
        }
    }

    boolean skip9 = FALSE;
    float a6 = yAFlipMask0 + a0;
    if (a6 > 0) {
        float x6 = x0;
        float y6 = y0 - (yNMask | 1);
        float z6 = z0;
        value += (a6 * a6) * (a6 * a6) * noise3_grad(seed,
                xrbp + (xNMask & S_PRIME_X), yrbp + (~yNMask & S_PRIME_Y), zrbp + (zNMask & S_PRIME_Z), x6, y6, z6);
    }
    else
    {
        float a7 = xAFlipMask0 + zAFlipMask0 + a0;
        if (a7 > 0) {
            float x7 = x0 - (xNMask | 1);
            float y7 = y0;
            float z7 = z0 - (zNMask | 1);
            value += (a7 * a7) * (a7 * a7) * noise3_grad(seed,
                    xrbp + (~xNMask & S_PRIME_X), yrbp + (yNMask & S_PRIME_Y), zrbp + (~zNMask & S_PRIME_Z), x7, y7, z7);
        }

        float a8 = yAFlipMask1 + a1;
        if (a8 > 0) {
            float x8 = x1;
            float y8 = (yNMask | 1) + y1;
            float z8 = z1;
            value += (a8 * a8) * (a8 * a8) * noise3_grad(seed2,
                    xrbp + S_PRIME_X, yrbp + (yNMask & (S_PRIME_Y << 1)), zrbp + S_PRIME_Z, x8, y8, z8);
            skip9 = TRUE;
        }
    }

    boolean skipD = FALSE;
    float aA = zAFlipMask0 + a0;
    if (aA > 0) {
        float xA = x0;
        float yA = y0;
        float zA = z0 - (zNMask | 1);
        value += (aA * aA) * (aA * aA) * noise3_grad(seed,
                xrbp + (xNMask & S_PRIME_X), yrbp + (yNMask & S_PRIME_Y), zrbp + (~zNMask & S_PRIME_Z), xA, yA, zA);
    }
    else
    {
        float aB = xAFlipMask0 + yAFlipMask0 + a0;
        if (aB > 0) {
            float xB = x0 - (xNMask | 1);
            float yB = y0 - (yNMask | 1);
            float zB = z0;
            value += (aB * aB) * (aB * aB) * noise3_grad(seed,
                    xrbp + (~xNMask & S_PRIME_X), yrbp + (~yNMask & S_PRIME_Y), zrbp + (zNMask & S_PRIME_Z), xB, yB, zB);
        }

        float aC = zAFlipMask1 + a1;
        if (aC > 0) {
            float xC = x1;
            float yC = y1;
            float zC = (zNMask | 1) + z1;
            value += (aC * aC) * (aC * aC) * noise3_grad(seed2,
                    xrbp + S_PRIME_X, yrbp + S_PRIME_Y, zrbp + (zNMask & (S_PRIME_Z << 1)), xC, yC, zC);
            skipD = TRUE;
        }
    }

    if (!skip5) {
        float a5 = yAFlipMask1 + zAFlipMask1 + a1;
        if (a5 > 0) {
            float x5 = x1;
            float y5 = (yNMask | 1) + y1;
            float z5 = (zNMask | 1) + z1;
            value += (a5 * a5) * (a5 * a5) * noise3_grad(seed2,
                    xrbp + S_PRIME_X, yrbp + (yNMask & (S_PRIME_Y << 1)), zrbp + (zNMask & (S_PRIME_Z << 1)), x5, y5, z5);
        }
    }

    if (!skip9) {
        float a9 = xAFlipMask1 + zAFlipMask1 + a1;
        if (a9 > 0) {
            float x9 = (xNMask | 1) + x1;
            float y9 = y1;
            float z9 = (zNMask | 1) + z1;
            value += (a9 * a9) * (a9 * a9) * noise3_grad(seed2,
                    xrbp + (xNMask & ((unsigned long)S_PRIME_X * 2)), yrbp + S_PRIME_Y, zrbp + (zNMask & (S_PRIME_Z << 1)), x9, y9, z9);
        }
    }

    if (!skipD) {
        float aD = xAFlipMask1 + yAFlipMask1 + a1;
        if (aD > 0) {
            float xD = (xNMask | 1) + x1;
            float yD = (yNMask | 1) + y1;
            float zD = z1;
            value += (aD * aD) * (aD * aD) * noise3_grad(seed2,
                    xrbp + (xNMask & (S_PRIME_X << 1)), yrbp + (yNMask & (S_PRIME_Y << 1)), zrbp + S_PRIME_Z, xD, yD, zD);
        }
    }

    return value;
}

/**
 * 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Y).
 * Recommended for 3D terrain and time-varied animations.
 * The Z coordinate should always be the "different" coordinate in whatever your use case is.
 * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, z, Y) or use noise3_XZBeforeY.
 * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, y, Z).
 * For a time varied animation, call noise3_ImproveXY(x, y, T).
 */
float OpenSimplex2S_noise3_ImproveXY(long seed, double x, double y, double z) {

    // Re-orient the cubic lattices without skewing, so Z points up the main lattice diagonal,
    // and the planes formed by XY are moved far out of alignment with the cube faces.
    // Orthonormal rotation. Not a skew transform.
    double xy = x + y;
    double s2 = xy * S_ROTATE3_ORTHOGONALIZER;
    double zz = z * S_ROOT3OVER3;
    double xr = x + s2 + zz;
    double yr = y + s2 + zz;
    double zr = xy * -S_ROOT3OVER3 + zz;

    // Evaluate both lattices to form a BCC lattice.
    return noise3_UnrotatedBase(seed, xr, yr, zr);
}

/**
 * 3D OpenSimplex2S/SuperSimplex noise, with better visual isotropy in (X, Z).
 * Recommended for 3D terrain and time-varied animations.
 * The Y coordinate should always be the "different" coordinate in whatever your use case is.
 * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, Y, z).
 * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, Z, y) or use noise3_ImproveXY.
 * For a time varied animation, call noise3_ImproveXZ(x, T, y) or use noise3_ImproveXY.
 */
float OpenSimplex2S_noise3_ImproveXZ(long seed, double x, double y, double z) {

    // Re-orient the cubic lattices without skewing, so Y points up the main lattice diagonal,
    // and the planes formed by XZ are moved far out of alignment with the cube faces.
    // Orthonormal rotation. Not a skew transform.
    double xz = x + z;
    double s2 = xz * -0.211324865405187;
    double yy = y * S_ROOT3OVER3;
    double xr = x + s2 + yy;
    double zr = z + s2 + yy;
    double yr = xz * -S_ROOT3OVER3 + yy;

    // Evaluate both lattices to form a BCC lattice.
    return noise3_UnrotatedBase(seed, xr, yr, zr);
}

/**
 * 3D OpenSimplex2S/SuperSimplex noise, fallback rotation option
 * Use noise3_ImproveXY or noise3_ImproveXZ instead, wherever appropriate.
 * They have less diagonal bias. This function's best use is as a fallback.
 */
float OpenSimplex2S_noise3_Fallback(long seed, double x, double y, double z) {

    // Re-orient the cubic lattices via rotation, to produce a familiar look.
    // Orthonormal rotation. Not a skew transform.
    double r = S_FALLBACK_ROTATE3 * (x + y + z);
    double xr = r - x, yr = r - y, zr = r - z;

    // Evaluate both lattices to form a BCC lattice.
    return noise3_UnrotatedBase(seed, xr, yr, zr);
}

/**
 * 4D SuperSimplex noise base.
 * Using ultra-simple 4x4x4x4 lookup partitioning.
 * This isn't as elegant or SIMD/GPU/etc. portable as other approaches,
 * but it competes performance-wise with optimized 2014 OpenSimplex.
 */
static float noise4_UnskewedBase(long seed, double xs, double ys, double zs, double ws) {

    // Get base points and offsets
    int xsb = fastFloor(xs), ysb = fastFloor(ys), zsb = fastFloor(zs), wsb = fastFloor(ws);
    float xsi = (float)(xs - xsb), ysi = (float)(ys - ysb), zsi = (float)(zs - zsb), wsi = (float)(ws - wsb);

    // Unskewed offsets
    float ssi = (xsi + ysi + zsi + wsi) * S_UNSKEW_4D;
    float xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;

    // Prime pre-multiplication for hash.
    long xsvp = xsb * S_PRIME_X, ysvp = ysb * S_PRIME_Y, zsvp = zsb * S_PRIME_Z, wsvp = wsb * S_PRIME_W;

    // Index into initial table.
    int index = ((fastFloor(xs * 4) & 3) << 0)
            | ((fastFloor(ys * 4) & 3) << 2)
            | ((fastFloor(zs * 4) & 3) << 4)
            | ((fastFloor(ws * 4) & 3) << 6);

    // Point contributions
    float value = 0;
    int secondaryIndexStartAndStop = LOOKUP_4D_A[index];
    int secondaryIndexStart = secondaryIndexStartAndStop & 0xFFFF;
    int secondaryIndexStop = secondaryIndexStartAndStop >> 16;
    for (int i = secondaryIndexStart; i < secondaryIndexStop; i++) {
        LatticeVertex4D c = LOOKUP_4D_B[i];
        float dx = xi + c.dx, dy = yi + c.dy, dz = zi + c.dz, dw = wi + c.dw;
        float a = (dx * dx + dy * dy) + (dz * dz + dw * dw);
        if (a < S_RSQUARED_4D) {
            a -= S_RSQUARED_4D;
            a *= a;
            value += a * a * noise4_grad(seed, xsvp + c.xsvp, ysvp + c.ysvp, zsvp + c.zsvp, wsvp + c.wsvp, dx, dy, dz, dw);
        }
    }
    return value;
}

/**
 * 4D SuperSimplex noise, with XYZ oriented like noise3_ImproveXY
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 * in a space where Z is vertical
 */
float OpenSimplex2S_noise4_ImproveXYZ_ImproveXY(long seed, double x, double y, double z, double w) {
    double xy = x + y;
    double s2 = xy * -0.21132486540518699998;
    double zz = z * 0.28867513459481294226;
    double ww = w * 1.118033988749894;
    double xr = x + (zz + ww + s2), yr = y + (zz + ww + s2);
    double zr = xy * -0.57735026918962599998 + (zz + ww);
    double wr = z * -0.866025403784439 + ww;

    return noise4_UnskewedBase(seed, xr, yr, zr, wr);
}

/**
 * 4D SuperSimplex noise, with XYZ oriented like noise3_ImproveXZ
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 * in a space where Y is vertical
 */
float OpenSimplex2S_noise4_ImproveXYZ_ImproveXZ(long seed, double x, double y, double z, double w) {
    double xz = x + z;
    double s2 = xz * -0.21132486540518699998;
    double yy = y * 0.28867513459481294226;
    double ww = w * 1.118033988749894;
    double xr = x + (yy + ww + s2), zr = z + (yy + ww + s2);
    double yr = xz * -0.57735026918962599998 + (yy + ww);
    double wr = y * -0.866025403784439 + ww;

    return noise4_UnskewedBase(seed, xr, yr, zr, wr);
}

/**
 * 4D SuperSimplex noise, with XYZ oriented like noise3_Fallback
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 * where there isn't a clear distinction between horizontal and vertical
 */
float OpenSimplex2S_noise4_ImproveXYZ(long seed, double x, double y, double z, double w) {
    double xyz = x + y + z;
    double ww = w * 1.118033988749894;
    double s2 = xyz * -0.16666666666666666 + ww;
    double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

    return noise4_UnskewedBase(seed, xs, ys, zs, ws);
}

/**
 * 4D SuperSimplex noise, with XY and ZW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
 * Recommended for noise(x, y, sin(time), cos(time)) trick.
 */
float OpenSimplex2S_noise4_ImproveXY_ImproveZW(long seed, double x, double y, double z, double w) {
    
    double s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718;
    double t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942;
    double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;
    
    return noise4_UnskewedBase(seed, xs, ys, zs, ws);
}

/**
 * 4D SuperSimplex noise, fallback lattice orientation.
 */
float OpenSimplex2S_noise4_Fallback(long seed, double x, double y, double z, double w) {

    // Get points for A4 lattice
    double s = S_SKEW_4D * (x + y + z + w);
    double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

    return noise4_UnskewedBase(seed, xs, ys, zs, ws);
}
