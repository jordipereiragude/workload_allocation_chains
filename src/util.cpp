/*
 * Copyright (c) 2023 Jordi Pereira, Marcus Ritt
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "util.hpp"

using namespace std;

options_counter verbosec;
State S;

string get_date_string(Timepoint t) {
#ifdef IS_MACOS
  char some_buffer[128];
  return string{some_buffer}; //macOS implementation of to_time_t is different to the gcc implementation (there are issues due to time unit and use of auto)
#else
  auto as_time_t = chrono::system_clock::to_time_t(t);
  struct tm tm;
  char some_buffer[128];
  if (::gmtime_r(&as_time_t, &tm))
    if (strftime(some_buffer, sizeof(some_buffer), "%F %T", &tm))
      return string{some_buffer};
  throw runtime_error("Failed to get current date as string");
#endif
}
