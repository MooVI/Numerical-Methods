/* 
 * File:   DESolvers.h
 * Author: qfeuille
 *
 * Created on December 21, 2011, 7:22 PM
 */

#ifndef DESOLVERS_H
#define	DESOLVERS_H



#include "FixedSteppers.h"
#include "VarSteppers.h"

#include "PRKMethods.h"
#include "RKMethods.h"
#include "RKNMethods.h"
#include "MultiSolvers.h"
#include "Numerovlike.h"

#include "Interrupts.h"

namespace NumMethod {
 inline int getnumSteps(double Stepsize, double begin, double end) {
        return ((int) ceil((end - begin) / Stepsize));
    }
}

#endif	/* DESOLVERS_H */

