#ifndef OMEGA_H_TEUCHOS_INCLUDES_HPP
#define OMEGA_H_TEUCHOS_INCLUDES_HPP

#include <Omega_h_macros.h>

OMEGA_H_SYSTEM_HEADER

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <Teuchos_Comm.hpp>
#include <Teuchos_Language.hpp>
#include <Teuchos_MathExpr.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Reader.hpp>
#include <Teuchos_TwoDArray.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_any.hpp>

#include <Teuchos_YamlParameterListHelpers.hpp>

#ifdef OMEGA_H_USE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#endif
