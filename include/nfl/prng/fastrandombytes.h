/*
 * File:   lattisigns512-20130329/fastrandombytes.h
 * Author: Gim Güneysu, Tobias Oder, Thomas Pöppelmann, Peter Schwabe
 * Public Domain
 */

#ifndef FASTRANDOMBYTES_H
#define FASTRANDOMBYTES_H

namespace nfl {
void fastrandombytes_seed(const unsigned char *s);
void fastrandombytes_reseed();
void fastrandombytes(unsigned char *r, unsigned long long rlen);
}

#endif
