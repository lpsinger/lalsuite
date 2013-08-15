## Check SWIG Octave module wrapping lal
## Author: Karl Wette, 2011, 2012

expected_exception = 0;

## check module load
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
disp("passed module load");

## check memory allocation
if lalcvar.swig_debug
  CheckMemoryLeaks();
  mem1 = new_Detector();
  mem2 = CreateCOMPLEX8Vector(5);
  mem3 = CreateREAL8Vector(3);
  mem4 = CreateREAL4TimeSeries("test", LIGOTimeGPS(0), 100, 0.1, lalcvar.lalDimensionlessUnit, 10);
  disp("*** below should be an error message from CheckMemoryLeaks() ***");
  try
    CheckMemoryLeaks();
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
  disp("*** above should be an error message from CheckMemoryLeaks() ***");
  clear mem1 mem2 mem3 mem4;
  CheckMemoryLeaks();
  disp("passed memory allocation");
else
  disp("skipped memory allocation");
endif

## check string conversions
strs = {"a"; "bc"; "def"};
sv = CreateStringVector(strs{:});
assert(sv.length == 3);
assert(all(strcmp(sv.data, strs)));
strs{1} = "ghijk";
sv.data{1} = strs{1};
strs{end+1} = "lmnopq";
sv = AppendString2Vector(sv, strs{4});
assert(sv.length == 4);
for i = 1:4
  assert(strcmp(sv.data{i}, strs{i}));
endfor
clear sv;
CheckMemoryLeaks();
disp("passed string conversions");

## check static vector/matrix conversions
lalcvar.lalswig_test_struct_vector{1} = lalcvar.lalswig_test_struct_const;
assert(lalcvar.lalswig_test_struct_vector{1}.n == lalcvar.lalswig_test_struct_const.n);
assert(lalcvar.lalswig_test_struct_vector{1}.i == lalcvar.lalswig_test_struct_const.i);
assert(lalcvar.lalswig_test_struct_vector{1}.f == lalcvar.lalswig_test_struct_const.f);
assert(strcmp(lalcvar.lalswig_test_struct_vector{1}.str, lalcvar.lalswig_test_struct_const.str));
assert(all(lalcvar.lalswig_test_struct_vector{1}.vec == lalcvar.lalswig_test_struct_const.vec));
lalcvar.lalswig_test_struct_matrix{1, 1} = lalcvar.lalswig_test_struct_const;
assert(lalcvar.lalswig_test_struct_matrix{1, 1}.n == lalcvar.lalswig_test_struct_const.n);
assert(lalcvar.lalswig_test_struct_matrix{1, 1}.i == lalcvar.lalswig_test_struct_const.i);
assert(lalcvar.lalswig_test_struct_matrix{1, 1}.f == lalcvar.lalswig_test_struct_const.f);
assert(strcmp(lalcvar.lalswig_test_struct_matrix{1, 1}.str, lalcvar.lalswig_test_struct_const.str));
assert(all(lalcvar.lalswig_test_struct_matrix{1, 1}.vec == lalcvar.lalswig_test_struct_const.vec));
sts = new_lalswig_test_struct();
assert(length(sts.vec) == 3);
assert(length(sts.evec) == 3);
assert(all(size(sts.mat) == [2, 3]));
sts.vec = [3; 2; 1];
assert(all(sts.vec == [3; 2; 1]));
sts.mat = [4, 5, 6; 9, 8, 7];
try
  sts.mat = [1.1, 2.3, 4.5; 6.5, 4.3, 2.1];
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
assert(all(all(sts.mat == [4, 5, 6; 9, 8, 7])));
for i = 1:3
  sts.evec(i) = 2*i + 3;
  assert(sts.evec(i) == (2*i + 3));
endfor
clear sts;
assert(!any(lalcvar.lalswig_test_enum_vector));
assert(!any(lalcvar.lalswig_test_enum_matrix(:)));
assert(length(lalcvar.lalswig_test_empty_INT4_vector) == 0);
assert(!any(lalcvar.lalswig_test_INT4_vector));
assert(!any(lalcvar.lalswig_test_INT4_matrix(:)));
assert(!any(lalcvar.lalswig_test_REAL8_vector));
assert(!any(lalcvar.lalswig_test_REAL8_matrix(:)));
assert(!any(lalcvar.lalswig_test_COMPLEX8_vector));
assert(!any(lalcvar.lalswig_test_COMPLEX8_matrix(:)));
lalcvar.lalswig_test_INT4_vector(1) = 10;
assert(lalcvar.lalswig_test_INT4_vector(1) == 10);
lalcvar.lalswig_test_INT4_matrix(1, 1) = 11;
assert(lalcvar.lalswig_test_INT4_matrix(1, 1) == 11);
lalcvar.lalswig_test_INT4_vector = lalcvar.lalswig_test_INT4_const_vector;
assert(all(lalcvar.lalswig_test_INT4_vector == [1; 2; 4]));
assert(lalcvar.lalswig_test_INT4_const_vector(3) == 4);
lalcvar.lalswig_test_INT4_matrix = lalcvar.lalswig_test_INT4_const_matrix;
assert(all(lalcvar.lalswig_test_INT4_matrix == [[1, 2, 4]; [2, 4, 8]]));
assert(lalcvar.lalswig_test_INT4_const_matrix(2, 3) == 8);
try
  lalcvar.lalswig_test_INT4_const_vector(20);
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
lalcvar.lalswig_test_REAL8_vector(1) = 3.4;
assert(lalcvar.lalswig_test_REAL8_vector(1) == 3.4);
lalcvar.lalswig_test_REAL8_matrix(1, 1) = 5.6;
assert(lalcvar.lalswig_test_REAL8_matrix(1, 1) == 5.6);
lalcvar.lalswig_test_COMPLEX8_vector(1) = complex(3.5, 4.75);
assert(lalcvar.lalswig_test_COMPLEX8_vector(1) == complex(3.5, 4.75));
lalcvar.lalswig_test_COMPLEX8_matrix(1, 1) = complex(5.5, 6.25);
assert(lalcvar.lalswig_test_COMPLEX8_matrix(1, 1) == complex(5.5, 6.25));
disp("passed static vector/matrix conversions");

## check dynamic vector/matrix conversions
function check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2)
  expected_exception = 0;
  iv.data = zeros(ivl, 1);
  rv.data = zeros(rvl, 1);
  cm.data = zeros(cms1, cms2);
  assert(ivl == 5);
  iv.data = [1; 3; 2; 4; 3];
  assert(all(iv.data == [1; 3; 2; 4; 3]));
  iv.data(4) = 7;
  assert(iv.data(4) == 7);
  assert(rvl == 5);
  rv.data = [1.2; 3.4; 2.6; 4.8; 3.5];
  assert(all(rv.data == [1.2; 3.4; 2.6; 4.8; 3.5]));
  rv.data(rvl) = 7.5;
  assert(rv.data(rvl) == 7.5);
  try
    rv.data(rvl + 1) = 99.9;
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
  try
    iv.data = rv.data;
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
  rv.data = iv.data;
  assert(all(rv.data == iv.data));
  assert(cms1 == 4);
  assert(cms2 == 6);
  for i = 1:cms1
    for j = 1:cms2
      cm.data(i, j) = complex(i / 4.0, j / 2.0);
    endfor
  endfor
  assert(cm.data(2, 3) == complex(0.5, 1.5));
  assert(cm.data(3, 2) == complex(0.75, 1.0));
  try
    iv.data(0) = cm.data(2, 3);
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
  try
    rv.data(0) = cm.data(3, 2);
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
endfunction
## check LAL vector and matrix datatypes
iv = CreateINT4Vector(5);
rv = CreateREAL8Vector(5);
cm = CreateCOMPLEX8VectorSequence(4, 6);
check_dynamic_vector_matrix(iv, iv.length, rv, rv.length,
                            cm, cm.length, cm.vectorLength);
clear iv rv cm;
rv0 = CreateREAL8Vector(0);
assert(rv0.length == 0);
assert(length(rv0.data) == 0);
clear rv0;
rv1 = CreateREAL8Vector(1);
rv1.data(1) = 1;
clear rv1;
CheckMemoryLeaks();
disp("passed dynamic vector/matrix conversions (LAL)");
## check GSL vectors and matrices
iv = new_gsl_vector_int(5);
rv = new_gsl_vector(5);
cm = new_gsl_matrix_complex_float(4, 6);
check_dynamic_vector_matrix(iv, iv.size, rv, rv.size,
                            cm, cm.size1, cm.size2);
clear iv rv cm;
rv1 = new_gsl_vector(1);
rv1.data(1) = 1;
clear rv1;
disp("passed dynamic vector/matrix conversions (GSL)");

## check dynamic array of pointers access
ap = lalswig_test_Create_arrayofptrs(3);
assert(ap.length == 3);
for i = 1:ap.length
  assert(ap.data{i}.length == 6);
  for j = 1:ap.data{i}.length
    assert(ap.data{i}.data(j) == 42*ap.length*(i-1) + (j-1));
  endfor
endfor
clear ap;
CheckMemoryLeaks();
disp("passed dynamic array of pointers access");

## check 'tm' struct conversions
gps = 989168284;
utc = [2011, 5, 11, 16, 57, 49, 4, 131, 0];
assert(all(GPSToUTC(gps) == utc));
assert(UTCToGPS(utc) == gps);
assert(UTCToGPS(utc(1:6)) == gps);
utc(7) = utc(8) = 0;
for i = [-1, 0, 1]
  utc(9) = i;
  assert(UTCToGPS(utc) == gps);
endfor
utcd = utc;
for i = 0:9
  utcd(3) = utc(3) + i;
  utcd = GPSToUTC(UTCToGPS(utcd));
  assert(utcd(7) == weekday(datenum(utcd(1:6))));
endfor
CheckMemoryLeaks();
disp("passed 'tm' struct conversions");

## check LIGOTimeGPS operations
t0 = new_LIGOTimeGPS();
assert(t0 == 0 && strcmp(swig_type(t0), "LIGOTimeGPS"));
t1 = new_LIGOTimeGPS(10.5);
t2 = new_LIGOTimeGPS(10, 500000000);
assert(t1 == t2 && strcmp(swig_type(t1), "LIGOTimeGPS"));
t3 = +t1;
t3 = -t2;
assert(t1 == t2 && t1 >= t2 && t2 >= t1);
assert(t1 + 3.5 == 14 && strcmp(swig_type(t1 + 3.5), "LIGOTimeGPS"));
assert(3.5 + t1 == 14 && strcmp(swig_type(3.5 + t1), "LIGOTimeGPS"));
t2 -= 5.5;
assert(t2 == 5 && strcmp(swig_type(t2), "LIGOTimeGPS"));
assert(t2 + 5.5 >= t1 && t2 + 3 != t2);
assert(t2 - 5 == t0 && strcmp(swig_type(t2 - 5), "LIGOTimeGPS"));
assert(t1 * 3 == 31.5 && strcmp(swig_type(t1 * 3), "LIGOTimeGPS"));
assert(3 * t1 == 31.5 && strcmp(swig_type(3 * t1), "LIGOTimeGPS"));
assert(t2 / 2.5 == 2 && strcmp(swig_type(t2 / 2.5), "LIGOTimeGPS"));
assert(21 / t1  == 2 && strcmp(swig_type(21 / t1), "LIGOTimeGPS"));
assert(t1 + t2 == 15.5 && strcmp(swig_type(t1 + t2), "LIGOTimeGPS"));
assert(t1 - t2 == 5.5 && strcmp(swig_type(t1 - t2), "LIGOTimeGPS"));
assert(t1 * t2 == 52.5 && strcmp(swig_type(t1 * t2), "LIGOTimeGPS"));
assert(t2 * t1 == 52.5 && strcmp(swig_type(t2 * t1), "LIGOTimeGPS"));
assert(t1 / t2 == 2.1 && strcmp(swig_type(t1 / t2), "LIGOTimeGPS"));
assert(t1 > t2 && t2 < t1 && t1 >= t2 && t2 <= t1)
assert(LIGOTimeGPS(333333333,333333333) == LIGOTimeGPS(1000000000) / 3)
assert(LIGOTimeGPS(666666666,666666667) == LIGOTimeGPS(2000000000) / 3)
t1 += 812345667.75;
assert(strcmp(t1.__str__(), "812345678.250000000"));
assert(new_LIGOTimeGPS(t1.__repr__()) == t1);
assert(t1.ns() == 812345678250000000);
t4struct = new_lalswig_test_gps;
t4struct.t = 1234.5;
assert(t4struct.t == 1234.5);
t5 = LIGOTimeGPS("1000");
assert(t5 == 1000);
if lalcvar.swig_version >= hex2dec("020010")
try
  t5 = LIGOTimeGPS("abc1000");
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
try
  t5 = LIGOTimeGPS("1000abc");
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
endif # lalcvar.swig_version >= hex2dec("020010")
assert(lalswig_test_noptrgps(LIGOTimeGPS(1234.5)) == lalswig_test_noptrgps(1234.5))
clear t0 t1 t2 t3 t4struct t5;
CheckMemoryLeaks();
disp("passed LIGOTimeGPS operations");

# check object parent tracking
a = lal.new_gsl_vector(3);
a.data = [1.1; 2.2; 3.3];
b = a.data;
assert(strcmp(typeinfo(b), "swiglal_oct_array_view_double"));
assert(all(b == [1.1; 2.2; 3.3]));
clear a;
assert(all(b == [1.1; 2.2; 3.3]));
ts = lal.CreateREAL8TimeSeries("test", lal.LIGOTimeGPS(0), 0, 0.1, lalcvar.lalDimensionlessUnit, 10);
ts.data.data = (1:10)';
for i = 1:7
  v = ts.data;
endfor
assert(v.data == (1:10)');
clear ts;
assert(v.data == (1:10)');
clear v;
CheckMemoryLeaks();
disp("passed object parent tracking");

## passed all tests!
disp("PASSED all tests");