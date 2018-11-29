/*
 * dyndist/tests/main.cpp
 *
 * Copyright (c) 2014,2018 Florian Pflug.
 *
 * This file is part of dyndist
 *
 * dyndist is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * dyndist is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with dyndist.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "dyndist/global.h"

#include <stdlib.h>

#include <stdexcept>

DYNDIST_NOWARN_PUSH
DYNDIST_NOWARN_SIGNS
#define BOOST_TEST_MODULE dyndist
#include <boost/test/unit_test.hpp>
DYNDIST_NOWARN_POP

