// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// tests for the BlockMask class
//
// here: BlockMask::operator|


#include <deal.II/fe/block_mask.h>

#include "../tests.h"



void
test()
{
  std::vector<bool> v1(12);
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1[i] = (i % 3 == 0);
  std::vector<bool> v2(12);
  for (unsigned int i = 0; i < v2.size(); ++i)
    v2[2] = (i % 4 == 0);

  BlockMask m1(v1);
  BlockMask m2(v2);
  BlockMask m = m1 | m2;

  // verify equality
  for (unsigned int i = 0; i < v1.size(); ++i)
    AssertThrow(m[i] == (v1[i] || v2[i]), ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
