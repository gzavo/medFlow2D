/*
 * Copyright (c) 2014, All Right Reserved, Gábor Závodszky, gabor@zavodszky.com

 * This source is subject to the
 * Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
 * Please see the License.txt file for more information.
 * All other rights reserved.

 * You should have received a copy of the license along with this
 * work in the License.txt file.
 * If not, see <http://creativecommons.org/licenses/by-nc-nd/4.0/>.
*/


#ifndef COMMON_H
#define COMMON_H

#include <stdbool.h>

typedef double real;

#ifndef M_PI
#	define M_PI 3.14159265358979323846
#endif

//#define PARALLEL
#define NUM_THREADS 8

#define my_free(x) do { free(x); x = NULL; } while (0)

#endif
