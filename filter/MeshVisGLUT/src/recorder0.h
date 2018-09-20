/*

Defines class recorder<Timer> for recording computing times as
measured by objects of class Timer.  See also recorder.h, which
defines another recorder class capable of also recording operation
counts.

*/

/*
 * Copyright (c) 1997 Rensselaer Polytechnic Institute
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Rensselaer Polytechnic Institute makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */

#ifndef RECORDER_H
#define RECORDER_H

#include <vector.h>

template <class Timer>
class recorder {
  vector<double> times;

public:
    void record(const Timer& t) {
      times.push_back(t.time());
  }
    void recordTime(const double& t) {
      times.push_back(t);
  }
  void report(ostream& o, int repeat_factor)
  {
    o << setiosflags(ios::fixed) << setprecision(3) << setw(12)
      << median(times)/repeat_factor;
      o << "\t";
  }
  void reset() {
    times.erase(times.begin(), times.end());
  }
  double median( vector<double> times) {
	  vector<double>::iterator midpoint = times.begin() + (times.end() - times.begin())/2;
	  nth_element(times.begin(), midpoint, times.end());
	  return *midpoint;

  }
};

#endif
