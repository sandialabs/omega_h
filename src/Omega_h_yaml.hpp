#ifndef OMEGA_H_YAML_HPP
#define OMEGA_H_YAML_HPP

#include <Omega_h_language.hpp>
#include <Omega_h_reader_tables.hpp>

namespace Omega_h {
namespace yaml {

Language build_language();
LanguagePtr ask_language();
ReaderTablesPtr ask_reader_tables();

}  // end namespace yaml
}  // end namespace Omega_h

#endif
