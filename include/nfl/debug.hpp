/* Copyright (C) 2015  Carlos Aguilar, Tancrède Lepoint, Adrien Guinet and Serge Guelton
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
#ifndef NFL_DEBUG_HPP
#define NFL_DEBUG_HPP

/***
 * Debug macros for NFL
 *
 * Assertions should be enabled for the development phase
 * Disable them after the initial development to tune performance
 * Disabling CHECK_STRICTMOD should improve performance automatically
 * Without ENFORCE_STRICTMOD arithmetic operations double the operand
 * size, be cautious with overflows in this case. As moduli are at least
 * two bits smaller than limbs there is some margin to increase
 * performance significantly.
 */

#include <cassert>

// ENFORCE_STRICTMOD in NTT this SHOULD NOT be deactivated
#define NTT_STRICTMOD

#ifdef CHECK_STRICTMOD
#define ASSERT_STRICTMOD(x) assert(x)
#else
#define ASSERT_STRICTMOD(x)
#endif

#endif
