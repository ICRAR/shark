//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
// Copyright by UWA (in the framework of the ICRAR)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

/**
 * @file
 *
 * The shark namespace description
 */

#ifndef INCLUDE_SHARK_H
#define INCLUDE_SHARK_H

/**
 * The shark namespace.
 *
 * The shark namespace contains all the code used internally by shark to run
 * its semi-analytical model. From a user's perspective only a handful of classes
 * are needed though:
 *
 * <ul>
 *  <li>The Options class, which encapsulates all user-provided options, and
 *  <li>The SharkRunner class, which takes an Options object, a number of threads,
 *      and runs shark until completion</li>
 * </ul>
 *
 * The @p shark binary uses these two classes to fire a shark instance. It
 * first parses the command-line parameters, constructs an Options objects and
 * determine the number of threads to use. It then creates a SharkRunner instance,
 * and finally invokes its @ref SharkRunner::run method.
 */
namespace shark {}

#endif // INCLUDE_SHARK_H