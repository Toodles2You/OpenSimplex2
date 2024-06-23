#ifndef __OPENSIMPLEX2_H
#define __OPENSIMPLEX2_H

#define PRIME_X 0x5205402B9270C86FL
#define PRIME_Y 0x598CD327003817B5L
#define PRIME_Z 0x5BCC226E9FA0BACBL
#define PRIME_W 0x56CC5227E58F554BL
#define HASH_MULTIPLIER 0x53A3F72DEEC546F5L
#define SEED_FLIP_3D -0x52D547B2E96ED629L
#define SEED_OFFSET_4D 0xE83DC3E0DA7164DL

#define ROOT2OVER2 0.7071067811865476
#define SKEW_2D 0.366025403784439
#define UNSKEW_2D -0.21132486540518713

#define ROOT3OVER3 0.577350269189626
#define FALLBACK_ROTATE_3D (2.0 / 3.0)
#define ROTATE_3D_ORTHOGONALIZER UNSKEW_2D

#define SKEW_4D -0.138196601125011f
#define UNSKEW_4D 0.309016994374947f
#define LATTICE_STEP_4D 0.2f

#define N_GRADS_2D_EXPONENT 7
#define N_GRADS_3D_EXPONENT 8
#define N_GRADS_4D_EXPONENT 9
#define N_GRADS_2D (1 << N_GRADS_2D_EXPONENT)
#define N_GRADS_3D (1 << N_GRADS_3D_EXPONENT)
#define N_GRADS_4D (1 << N_GRADS_4D_EXPONENT)

#define NORMALIZER_2D 0.01001634121365712
#define NORMALIZER_3D 0.07969837668935331
#define NORMALIZER_4D 0.0220065933241897

#define RSQUARED_2D 0.5f
#define RSQUARED_3D 0.6f
#define RSQUARED_4D 0.6f

#define GRADIENTS_2D_SRC "gradients/fast/GRADIENTS_2D.h"
#define GRADIENTS_3D_SRC "gradients/fast/GRADIENTS_3D.h"
#define GRADIENTS_4D_SRC "gradients/fast/GRADIENTS_4D.h"

#endif /* !__OPENSIMPLEX2_H */