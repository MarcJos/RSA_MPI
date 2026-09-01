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

// #include <exanb/core/grid_fields.h>
#ifndef XNB_HAS_GRID_FIELDS_DEFINTIONS
#error Cannot be included outside of exanb/core/grid_fields.h
#endif

#include <cstdint>

// position (rx, ry, rz) and id are already declared by exanb/core/grid_fields.h;
// rsa_data_storage's remaining per-sphere attributes become their own fields here.
XNB_DECLARE_FIELD(double, radius, "sphere radius");
XNB_DECLARE_FIELD(int32_t, phase, "sphere phase");
XNB_DECLARE_FIELD(uint64_t, priority, "draw priority, used to resolve placement conflicts");
XNB_DECLARE_FIELD(int32_t, confirmed,
                  "multi-pass candidate resolution: 1 once this candidate is proven to have no "
                  "surviving lower-priority neighbor (a local minimum, definitely accepted), 0 while "
                  "still undecided - see exanb_naive::resolve_candidates_pass");

namespace rsa_mpi {
// note: deliberately not "using namespace ::exanb" here - a using-directive
// inside "namespace rsa_mpi" leaks into every later reopening of that same
// namespace within the translation unit (e.g. any operator .cpp that ends up
// including this file transitively through exanb/core/grid_fields.h),
// making onika::scg's INPUT/OUTPUT/REQUIRED ambiguous with exanb's re-exported
// copies of the same symbols there.

// rx, ry and rz are added implicitly
using RSAFieldSet = ::exanb::FieldSet<::exanb::field::_id, ::exanb::field::_radius, ::exanb::field::_phase,
                                      ::exanb::field::_priority, ::exanb::field::_confirmed>;

static inline constexpr ::exanb::FieldSets<RSAFieldSet> available_field_sets_v = {};
}  // namespace rsa_mpi

#define HAS_POSITION_BACKUP_FIELDS false
#define PositionBackupFieldX ::exanb::unused_field_id_v
#define PositionBackupFieldY ::exanb::unused_field_id_v
#define PositionBackupFieldZ ::exanb::unused_field_id_v

#define XNB_AVAILABLE_FIELD_SETS ::rsa_mpi::available_field_sets_v
