/* Copyright (c) 2018 Julio Daniel Machado Silva */

/* Licensed under the Apache License, Version 2.0 (the "License"); */
/* you may not use this file except in compliance with the License. */
/* You may obtain a copy of the License at */

/* http://www.apache.org/licenses/LICENSE-2.0 */

/* Unless required by applicable law or agreed to in writing, software */
/* distributed under the License is distributed on an "AS IS" BASIS, */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/* See the License for the specific language governing permissions and */
/* limitations under the License. */

#include "nspline.hpp"

#include <iostream>
#include <vector>


namespace nspline {

double abs( double x )
{
    return x >= 0 ? x : -x;
}

bool test()
{
    // Input points { x_1, x_2, ..., x_n }
    std::vector<double> x = { 0.0, 0.5, 1.0 };

    // Input function samples { f(x_1), f(x_2), ..., f(x_n) }
    std::vector<double> y = { 1.0, 1.5, -1.0 };

    // Interpolant:
    jdms::NSpline f(x, y);

    bool success = true;
    double tolerance = 1E-5;

    std::cout << std::endl;

    for ( int i = 0; i < x.size(); ++i )
    {
        if ( abs(f( x[i] ) - y[i]) > tolerance )
        {
            std::cerr << "Error ---> `f(" << x[i] << ") = " << f(x[i]) 
                      << "` differs more than tolerance `" << tolerance 
                      << "` from input data `" << y[i] << "`" << std::endl;
            std::cout << std::endl << std::flush;

            return !success;
        }
    }

    std::cout << "All tests passed" << std::endl;
    std::cout << std::endl << std::flush;

    return success;
}

} // namespace nspline

int main()
{
    if ( !nspline::test() )
    {
        return 1;
    }

    return 0;
}
