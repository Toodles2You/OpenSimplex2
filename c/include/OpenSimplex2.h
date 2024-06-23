/**
 * K.jpg's OpenSimplex 2, faster variant
 */

#ifndef _OPENSIMPLEX2_H
#define _OPENSIMPLEX2_H

/**
 * 2D OpenSimplex2 noise.
 */
float OpenSimplex2_noise2(long seed, double x, double y);

/**
 * 2D Simplex noise, with Y pointing down the main diagonal.
 * Might be better for a 2D sandbox style game, where Y is vertical.
 * Probably slightly less optimal for heightmaps or continent maps,
 * unless your map is centered around an equator. It's a subtle
 * difference, but the option is here to make it an easy choice.
 */
float OpenSimplex2_noise2_ImproveX(long seed, double x, double y);

/**
 * 3D OpenSimplex2 noise, with better visual isotropy in (X, Y).
 * Recommended for 3D terrain and time-varied animations.
 * The Z coordinate should always be the "different" coordinate in whatever your use case is.
 * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, z, Y) or use noise3_XZBeforeY.
 * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, y, Z).
 * For a time varied animation, call noise3_ImproveXY(x, y, T).
 */
float OpenSimplex2_noise3_ImproveXY(long seed, double x, double y, double z);

/**
 * 3D OpenSimplex2 noise, with better visual isotropy in (X, Z).
 * Recommended for 3D terrain and time-varied animations.
 * The Y coordinate should always be the "different" coordinate in whatever your use case is.
 * If Y is vertical in world coordinates, call noise3_ImproveXZ(x, Y, z).
 * If Z is vertical in world coordinates, call noise3_ImproveXZ(x, Z, y) or use noise3_ImproveXY.
 * For a time varied animation, call noise3_ImproveXZ(x, T, y) or use noise3_ImproveXY.
 */
float OpenSimplex2_noise3_ImproveXZ(long seed, double x, double y, double z);

/**
 * 3D OpenSimplex2 noise, fallback rotation option
 * Use noise3_ImproveXY or noise3_ImproveXZ instead, wherever appropriate.
 * They have less diagonal bias. This function's best use is as a fallback.
 */
float OpenSimplex2_noise3_Fallback(long seed, double x, double y, double z);

/**
 * 4D OpenSimplex2 noise, with XYZ oriented like noise3_ImproveXY
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 * in a space where Z is vertical
 */
float OpenSimplex2_noise4_ImproveXYZ_ImproveXY(long seed, double x, double y, double z, double w);

/**
 * 4D OpenSimplex2 noise, with XYZ oriented like noise3_ImproveXZ
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 * in a space where Y is vertical
 */
float OpenSimplex2_noise4_ImproveXYZ_ImproveXZ(long seed, double x, double y, double z, double w);

/**
 * 4D OpenSimplex2 noise, with XYZ oriented like noise3_Fallback
 * and W for an extra degree of freedom. W repeats eventually.
 * Recommended for time-varied animations which texture a 3D object (W=time)
 * where there isn't a clear distinction between horizontal and vertical
 */
float OpenSimplex2_noise4_ImproveXYZ(long seed, double x, double y, double z, double w);

/**
 * 4D OpenSimplex2 noise, with XY and ZW forming orthogonal triangular-based planes.
 * Recommended for 3D terrain, where X and Y (or Z and W) are horizontal.
 * Recommended for noise(x, y, sin(time), cos(time)) trick.
 */
float OpenSimplex2_noise4_ImproveXY_ImproveZW(long seed, double x, double y, double z, double w);

/**
 * 4D OpenSimplex2 noise, fallback lattice orientation.
 */
float OpenSimplex2_noise4_Fallback(long seed, double x, double y, double z, double w);

#endif /* !_OPENSIMPLEX2_H */
