/*
   Licensed to the Apache Software Foundation (ASF) under one
   or more contributor license agreements.  See the NOTICE file
   distributed with this work for additional information
   regarding copyright ownership.  The ASF licenses this file
   to you under the Apache License, Version 2.0 (the
   "License"); you may not use this file except in compliance
   with the License.  You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
 */
#pragma once

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <yaml-cpp/yaml.h>
#include <memory>

#include <onika/type_utils.h>
#include <onika/log.h>

namespace rsa_mpi
{
  using namespace onika;
  using namespace scg;
  /*
     Internal template utilities
   */
  namespace details
  {
    template< template<int> typename _OperatorTemplate, int... DIM >
      struct MakeRSAMPIOperatorHelper
      {
        static inline std::shared_ptr<onika::scg::OperatorNode> make_operator( const YAML::Node& node, const onika::scg::OperatorNodeFlavor& flavor )
        {
          return make_compatible_operator < _OperatorTemplate<DIM>...> (node,flavor);
        }
      };

    template< template<int> class _OperatorTemplate, int... DIM>
      struct make_rsa_mpi_operator_t
      {
        static inline onika::scg::OperatorNodeCreateFunction make_factory(const std::string& opname)
        {
          onika::scg::OperatorNodeCreateFunction factory = [] (const YAML::Node& node, const onika::scg::OperatorNodeFlavor& flavor) -> std::shared_ptr<onika::scg::OperatorNode>
          {
            std::shared_ptr<onika::scg::OperatorNode> op = MakeRSAMPIOperatorHelper< _OperatorTemplate, DIM... >::make_operator(node,flavor);
            return op;        
          };

          return factory;
        }
      };

  } // temporary close exanb namespace
}

namespace onika
{
  namespace scg
  {
    template< template<int> class _OperatorTemplate, int... DIM>
      struct OperatorNodeFactoryGenerator< rsa_mpi::details::make_rsa_mpi_operator_t<_OperatorTemplate, DIM...> >
      {
        static inline OperatorNodeCreateFunction make_factory(const std::string& opname)
        {
          return rsa_mpi::details::make_rsa_mpi_operator_t<_OperatorTemplate, DIM...>::make_factory(opname) ;
        }
      };
  }
}

namespace rsa_mpi
{
  template< template<int> class _OperatorTemplate > 
    static inline constexpr 
    onika::scg::OperatorNodeFactoryGenerator< details::make_rsa_mpi_operator_t< _OperatorTemplate, 2, 3, 4, 5, 6, 7, 8, 9, 10 > > make_rsa_mpi_operator = {};
}

