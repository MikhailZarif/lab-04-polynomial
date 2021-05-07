// Copyright 2020 Your Name <your_email>

#include <gtest/gtest.h>

#include <polynomial.hpp>

TEST(Polynom, Init) {
  int* vec = new int[4];
  for (int i = 0; i < 4; i++) {
    vec[i] = i * i;
  }
  Polynom<int> a(4, vec);
  EXPECT_EQ(a[0], 0);
  EXPECT_EQ(a[1], 1);
  EXPECT_EQ(a[2], 4);
  EXPECT_EQ(a[3], 9);

  a * 5;
  EXPECT_EQ(a[0], 0);
  EXPECT_EQ(a[1], 5);
  EXPECT_EQ(a[2], 20);
  EXPECT_EQ(a[3], 45);
  delete[] vec;
}

TEST(Polynom, CopyAssign) {
  auto* vec = new double[4];
  for (int i = 0; i < 4; i++) {
    vec[i] = i + 0.1 * i;
  }
  Polynom<double> a(4, vec);
  Polynom<double> cp(a);
  EXPECT_EQ(cp, a);
  delete[] vec;

  a[0] = 0.25;
  a[2] = 1.5;
  EXPECT_NE(cp, a);

  cp = a;
  EXPECT_EQ(cp, a);
}

TEST(Polynom, Degree) {
  auto* vec = new int[5];
  for (int i = 0; i < 5; i++) {
    vec[i] = i;
  }
  Polynom<int> a(5, vec);
  size_t res = a.Degree();
  EXPECT_EQ(res, 4);
  delete[] vec;

  a[4] = 0;
  a[3] = 0;
  res = a.Degree();
  EXPECT_EQ(res, 2);
}

TEST(Polynom, Count_x) {
  auto* vec = new int[4];
  for (int i = 0; i < 4; i++) {
    vec[i] = i + 1;
  }
  // 4x^3 + 3x^2 + 2x + 1
  Polynom<int> a(4, vec);
  delete[] vec;
  int result = a.Func(0);
  EXPECT_EQ(result, 1);

  a[3] = -1;
  a[0] = 0;
  // -x^3 + 3x^2 +2x
  result = a.Func(10);
  EXPECT_EQ(result, -680);
}

TEST(Polynom, Add) {
  auto* vec1 = new int[4];
  for (int i = 0; i < 4; i++) {
    vec1[i] = i + 1;
  }
  auto* vec2 = new int[4];
  for (int i = 0; i < 4; i++) {
    vec2[i] = i + 1;
  }
  vec2[3] = -4;
  vec1[1] = -2;
  // (4x^3 + 3x^2 - 2x + 1) + (-4x^3 + 3x^2 + 2x + 1)
  Polynom<int> a(4, vec1);
  Polynom<int> b(4, vec2);
  delete[] vec1;
  delete[] vec2;
  // 6x^2 + 2
  Polynom<int> c = a + b;
  size_t deg = c.Degree();
  EXPECT_EQ(deg, 2);
  EXPECT_EQ(c[0], 2);
  EXPECT_EQ(c[1], 0);
  EXPECT_EQ(c[2], 6);

  c += a;
  // 4x^3 + 9x^2 - 2x + 3
  deg = c.Degree();
  EXPECT_EQ(deg, 3);
  EXPECT_EQ(c[0], 3);
  EXPECT_EQ(c[1], -2);
  EXPECT_EQ(c[2], 9);
  EXPECT_EQ(c[3], 4);

  Polynom<int> d = c;
  d * -1;
  c += d;
  deg = c.Degree();
  EXPECT_EQ(deg, 0);
}

TEST(Polynom, Sub) {
  auto* vec1 = new int[4];
  for (int i = 0; i < 4; i++) {
    vec1[i] = i + 1;
  }
  auto* vec2 = new int[4];
  for (int i = 0; i < 4; i++) {
    vec2[i] = i + 1;
  }
  // (4x^3 + 3x^2 + 2x + 1) - (4x^3 + 3x^2 + 2x + 1)
  Polynom<int> a(4, vec1);
  Polynom<int> b(4, vec2);
  delete[] vec1;
  delete[] vec2;
  Polynom<int> c = a - b;
  size_t deg = c.Degree();
  EXPECT_EQ(deg, 0);

  b[0] = -1;
  b[3] = 2;
  // 2x^3 + 2
  c = a - b;
  deg = c.Degree();
  EXPECT_EQ(deg, 3);
  EXPECT_EQ(c[0], 2);
  EXPECT_EQ(c[1], 0);
  EXPECT_EQ(c[2], 0);
  EXPECT_EQ(c[3], 2);

  c -= a;
  // -2x^3 - 3x^2 - 2x + 1
  deg = c.Degree();
  EXPECT_EQ(deg, 3);
  EXPECT_EQ(c[0], 1);
  EXPECT_EQ(c[1], -2);
  EXPECT_EQ(c[2], -3);
  EXPECT_EQ(c[3], -2);
}

TEST(Polynom, Mult) {
  auto* vec1 = new int[3];
  for (int i = 0; i < 3; i++) {
    vec1[i] = i + 1;
  }
  auto* vec2 = new int[2];
  for (int i = 0; i < 2; i++) {
    vec2[i] = i + 1;
  }
  Polynom<int> a(3, vec1);  // 3x^2 + 2x + 1
  Polynom<int> b(2, vec2);  // 2x + 1
  Polynom<int> c = a * b;   // 6x^3 + 7x^2 + 4x + 1
  delete[] vec1;
  delete[] vec2;
  size_t deg = c.Degree();
  EXPECT_EQ(deg, 3);
  EXPECT_EQ(c[0], 1);
  EXPECT_EQ(c[1], 4);
  EXPECT_EQ(c[2], 7);
  EXPECT_EQ(c[3], 6);
  a *= b;
  EXPECT_EQ(a, c);

  b[0] = 0;
  b[1] = 0;
  Polynom<int> d = b * c;
  deg = d.Degree();
  EXPECT_EQ(deg, 0);
  a *= b;
  deg = a.Degree();
  EXPECT_EQ(deg, 0);
  EXPECT_EQ(a, d);
}

TEST(Polynom, Div) {
  auto* vec1 = new int[3];
  vec1[2] = 1;
  vec1[1] = -2;
  vec1[0] = 1;
  auto* vec2 = new int[2];
  vec2[1] = 1;
  vec2[0] = -1;
  Polynom<int> a(3, vec1);    // x^2 - 2x + 1
  Polynom<int> b(2, vec2);    // x - 1
  Polynom<double> c = a / b;  // x - 1 обычное деление, на выходе только полином
  delete[] vec1;
  delete[] vec2;
  size_t deg = c.Degree();
  EXPECT_EQ(deg, 1);
  EXPECT_NEAR(c[0], -1, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(c[1], 1, std::numeric_limits<double>::epsilon());

  b[1] = 2;
  b[0] = 1;  // 2x + 1
  a[2] = 3;
  a[1] = 0;
  a[0] = -1;  // 3x^2 - 1
  c = a / b;  // 1.5x - 0.75 (остаток -0.25)
  deg = c.Degree();
  EXPECT_EQ(deg, 1);
  EXPECT_NEAR(c[0], -0.75, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(c[1], 1.5, std::numeric_limits<double>::epsilon());

  auto* v = new double[3];
  v[0] = -1;
  v[1] = 0;
  v[2] = 3;
  auto* m = new double[2];
  m[0] = 1;
  m[1] = 2;
  Polynom<double> db(2, m);
  Polynom<double> da(3, v);
  delete[] m;
  delete[] v;
  da /= db;  //Операция /= принимает и выводит только значения типа double,
             //чтобы не терять данные при смене типов
  deg = da.Degree();
  EXPECT_EQ(deg, 1);
  EXPECT_EQ(da, c);

  db[0] = 0;
  db[1] = 0;  //При делении на ноль вернём значение делимого, чтобы пользователь
              //не терял данные
  da /= db;
  deg = da.Degree();
  EXPECT_EQ(deg, 1);
  EXPECT_EQ(da, c);
}

TEST(Polynom, Mod) {
  auto* vec1 = new double[3];
  vec1[2] = 1;
  vec1[1] = -2;
  vec1[0] = 1;
  auto* vec2 = new double[2];
  vec2[1] = 1;
  vec2[0] = -1;
  Polynom<double> a(3, vec1);  // x^2 - 2x + 1
  Polynom<double> b(2, vec2);  // x - 1
  Polynom<double> c = a % b;   // остаток 0
  delete[] vec1;
  delete[] vec2;
  size_t deg = c.Degree();
  EXPECT_EQ(deg, 0);
  Polynom<double> d = a;
  a %= b;
  EXPECT_EQ(a, c);

  b[1] = 2;
  b[0] = 1;  // 2x + 1
  d[2] = 3;
  d[1] = 0;
  d[0] = -1;  // 3x^2 - 1
  c = d % b;  // остаток -0.25
  deg = c.Degree();
  EXPECT_EQ(deg, 0);
  EXPECT_NEAR(c[0], -0.25, std::numeric_limits<double>::epsilon());
  d %= b;
  EXPECT_EQ(c, d);

  a = c;
  a[0] = 0;
  c = d % a;  // Деление на 0
  EXPECT_EQ(c, d);
  d %= a;
  EXPECT_EQ(c, d);
  EXPECT_NEAR(c[0], -0.25, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(d[0], -0.25, std::numeric_limits<double>::epsilon());
}
