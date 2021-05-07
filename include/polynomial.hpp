// Copyright 2018 Your Name <your_email>
#include <cmath>
#include <limits>
#ifndef INCLUDE_POLYNOMIAL_HPP_
#define INCLUDE_POLYNOMIAL_HPP_
template <class T>
class Polynom {
  //Вектор коэффициентов
  //Для удобства deg - это не степень полинома а количество слагаемых в нём
 private:
  T* data;
  size_t deg;

 public:
  //Конструктор по умолчанию
  Polynom() {
    deg = 0;
    data = nullptr;
  }

  //Коструктор копирования
  Polynom(Polynom<T>& c) {
    deg = c.deg;
    if (deg != 0) {
      data = new T[deg];
      for (size_t i = 0; i < deg; i++) {
        data[i] = c[i];
      }
    } else {
      data = nullptr;
    }
  }

  //Пользовательский конструктор, где d - это количество слагаемых
  Polynom(const size_t d, const T* k) {
    if ((k == nullptr) || (d == 0)) {
      deg = 0;
      data = nullptr;
    } else {
      deg = d;
      data = new T[deg];
      for (size_t i = 0; i < deg; i++) {
        data[i] = k[i];
      }
    }
  }

  //Деструктор
  ~Polynom() { delete[] data; }

  //Оператор =
  Polynom& operator=(const Polynom& a) {
    if (&a != this) {
      delete[] data;
      deg = a.deg;
      if (deg != 0) {
        data = new T[a.deg];
        for (size_t i = 0; i < deg; i++) {
          data[i] = a.data[i];
        }
      } else {
        data = nullptr;
      }
    }
    return *this;
  }

  //Оператор сравнения ==
  bool operator==(const Polynom& a) const {
    bool flag = true;
    if (deg == a.deg) {
      if (std::is_floating_point<T>::value) {
        for (size_t i = 0; i < deg; i++) {
          if (abs(data[i] - a.data[i]) > std::numeric_limits<T>::epsilon()) {
            flag = false;
            break;
          }
        }
      } else {
        for (size_t i = 0; i < deg; i++) {
          if (data[i] != a.data[i]) {
            flag = false;
            break;
          }
        }
      }
    } else {
      flag = false;
    }
    return flag;
  }

  //Оператор сравнения !=
  bool operator!=(const Polynom& a) const {
    bool flag = false;
    if (deg == a.deg) {
      if (std::is_floating_point<T>::value) {
        for (size_t i = 0; i < deg; i++) {
          if (abs(data[i] - a.data[i]) > std::numeric_limits<T>::epsilon()) {
            flag = true;
            break;
          }
        }
      } else {
        for (size_t i = 0; i < deg; i++) {
          if (data[i] != a.data[i]) {
            flag = true;
            break;
          }
        }
      }
    } else {
      flag = true;
    }
    return flag;
  }

  //Получение и установка коэффициента при одночлене степени i
  T operator[](const size_t i) const { return data[i]; }
  T& operator[](const size_t i) { return (data[i]); }

  //Степень полинома
  [[nodiscard]] size_t Degree() const {
    size_t degree = 0;
    if (data != nullptr) {
      for (size_t i = deg - 1; i >= 0; i--) {
        if (data[i] != 0) {
          degree = i;
          break;
        }
        if (i == 0) {
          break;
        }
      }
    }
    return degree;
  }

  //Сложение +=
  Polynom<T>& operator+=(const Polynom<T>& b) {
    size_t newdeg;
    size_t min;
    bool flag = true;
    if (Degree() < b.Degree()) {
      newdeg = b.Degree() + 1;
      min = Degree() + 1;
    } else {
      newdeg = Degree() + 1;
      min = b.Degree() + 1;
      flag = false;
    }
    if (Degree() == b.Degree()) {
      size_t i = Degree();
      while ((data[i] + b[i] == 0) && (i >= 0) && (i <= Degree())) {
        if (i != 0) {
          newdeg--;
          min--;
          i--;
        } else {
          newdeg--;
          min--;
          break;
        }
      }
    }
    T* res;
    if (newdeg == 0) {
      res = nullptr;
    } else {
      res = new T[newdeg];
      for (size_t i = 0; i < min; i++) {
        res[i] = data[i] + b[i];
      }
      if (flag) {
        for (size_t j = min; j < newdeg; j++) {
          res[j] = b[j];
        }
      } else {
        for (size_t j = min; j < newdeg; j++) {
          res[j] = data[j];
        }
      }
    }
    delete[] data;
    deg = newdeg;
    if (newdeg != 0) {
      data = new T[newdeg];
      for (size_t i = 0; i < deg; i++) {
        data[i] = res[i];
      }
    } else {
      data = nullptr;
    }
    delete[] res;
    return *this;
  }

  //Вычитание -=
  Polynom<T>& operator-=(const Polynom<T>& b) {
    size_t newdeg;
    size_t min;
    bool flag = true;
    if (Degree() < b.Degree()) {
      newdeg = b.Degree() + 1;
      min = Degree() + 1;
    } else {
      newdeg = Degree() + 1;
      min = b.Degree() + 1;
      flag = false;
    }
    if (Degree() == b.Degree()) {
      size_t i = Degree();
      while ((data[i] - b[i] == 0) && (i >= 0) && (i <= Degree())) {
        if (i != 0) {
          newdeg--;
          min--;
          i--;
        } else {
          newdeg--;
          min--;
          break;
        }
      }
    }
    T* res;
    if (newdeg == 0) {
      res = nullptr;
    } else {
      res = new T[newdeg];
      for (size_t i = 0; i < min; i++) {
        res[i] = data[i] - b[i];
      }
      if (flag) {
        for (size_t j = min; j < newdeg; j++) {
          res[j] = -b[j];
        }
      } else {
        for (size_t j = min; j < newdeg; j++) {
          res[j] = data[j];
        }
      }
    }
    delete[] data;
    deg = newdeg;
    if (newdeg != 0) {
      data = new T[newdeg];
      for (size_t i = 0; i < deg; i++) {
        data[i] = res[i];
      }
    } else {
      data = nullptr;
    }
    delete[] res;
    return *this;
  }

  //Умножение на константу типа Т
  Polynom<T>& operator*(const T& d) {
    for (size_t i = 0; i < deg; i++) {
      data[i] *= d;
    }
    return *this;
  }

  //Умножение на полином *=
  Polynom<T>& operator*=(const Polynom<T>& a) {
    size_t sa = a.Degree();
    size_t sb = Degree();
    if (((a[0] == 0) && (sa == 0)) || ((data[0] == 0) && (sb == 0))) {
      delete[] data;
      data = nullptr;
      deg = 0;
    } else {
      size_t newdeg = sa + sb + 1;
      T* res = new T[newdeg];
      for (size_t k = 0; k < newdeg; k++) {
        res[k] = 0;
      }
      for (size_t i = 0; i < sa + 1; i++) {
        for (size_t j = 0; j < sb + 1; j++) {
          res[i + j] += a[i] * data[j];
        }
      }
      delete[] data;
      deg = newdeg;
      data = new T[deg];
      for (size_t i = 0; i < deg; i++) {
        data[i] = res[i];
      }
      delete[] res;
    }
    return *this;
  }

  //Деление на полином /= (принимает и выдаёт только double)
  Polynom<double>& operator/=(const Polynom<double>& b) {
    if (((data == nullptr) || ((Degree() == 0) && (data[0] == 0))) ||
        (((b.Degree() == 0) && (b[0] == 0)) || (b.data == nullptr))) {
      return *this;
    } else {
      auto* copyofa = new double[Degree() + 1];
      for (size_t i = 0; i < Degree() + 1; i++) {
        copyofa[i] = data[i];
      }
      Polynom<double> acp(Degree() + 1, copyofa);
      auto* copyofb = new double[b.Degree() + 1];
      for (size_t i = 0; i < b.Degree() + 1; i++) {
        copyofb[i] = b[i];
      }
      Polynom<double> bcp(b.Degree() + 1, copyofb);
      if (Degree() < b.Degree()) {
        delete[] data;
        deg = 1;
        data = new double[1];
        data[0] = 0;
        return *this;
      } else {
        double* res;
        res = new double[Degree() - b.Degree() + 1];
        double* min;
        min = new double[acp.Degree() - bcp.Degree() + 1];
        while (acp.Degree() >= bcp.Degree()) {
          size_t delimoe = acp.Degree();
          size_t delitel = bcp.Degree();
          res[delimoe - delitel] = acp[acp.Degree()] / bcp[bcp.Degree()];
          for (size_t i = 0; i < acp.Degree() - bcp.Degree() + 1; i++) {
            min[i] = 0;
          }
          min[delimoe - delitel] = res[delimoe - delitel];
          Polynom<double> m(delimoe - delitel + 1, min);
          m *= bcp;
          acp -= m;
        }
        deg = Degree() - b.Degree() + 1;
        delete[] data;
        data = new double[deg];
        for (size_t i = 0; i < deg; i++) {
          data[i] = res[i];
        }
        delete[] res;
        delete[] copyofa;
        delete[] copyofb;
        delete[] min;
        return *this;
      }
    }
  }

  //Остаток от деления на полином %=
  Polynom<double>& operator%=(const Polynom<double>& a) {
    if ((a.data == nullptr) || ((a.Degree() == 0) && (a[0] == 0))) {
      return *this;
    } else {
      auto* n = new double[deg];
      for (size_t i = 0; i < deg; i++) {
        n[i] = data[i];
      }
      Polynom<double> ch(deg, n);
      Polynom<double> cp(deg, n);
      ch /= a;
      ch *= a;
      cp -= ch;
      delete[] data;
      if (cp.data != nullptr) {
        deg = cp.Degree() + 1;
        data = new double[deg];
        for (size_t i = 0; i < deg; i++) {
          data[i] = cp[i];
        }
      } else {
        deg = 0;
        data = nullptr;
      }
      delete[] n;
      return *this;
    }
  }

  //Вычисление значения полинома при значении х
  T Func(const T x) {
    T arg = 1;
    T result = 0;
    for (size_t i = 0; i < deg; i++) {
      if (i != 0) {
        arg *= x;
      }
      result += data[i] * arg;
    }
    return result;
  }
};

//Сложение +
template <class T>
Polynom<T> operator+(const Polynom<T>& a, const Polynom<T>& b) {
  size_t newdeg;
  size_t min;
  bool flag = true;
  if (a.Degree() < b.Degree()) {
    newdeg = b.Degree() + 1;
    min = a.Degree() + 1;
  } else {
    newdeg = a.Degree() + 1;
    min = b.Degree() + 1;
    flag = false;
  }
  if (a.Degree() == b.Degree()) {
    size_t i = a.Degree();
    while ((a[i] + b[i] == 0) && (i >= 0) && (i <= a.Degree())) {
      if (i != 0) {
        newdeg--;
        min--;
        i--;
      } else {
        newdeg--;
        min--;
        break;
      }
    }
  }
  T* res;
  if (newdeg == 0) {
    res = nullptr;
  } else {
    res = new T[newdeg];
    for (size_t i = 0; i < min; i++) {
      res[i] = a[i] + b[i];
    }
    if (flag) {
      for (size_t j = min; j < newdeg; j++) {
        res[j] = b[j];
      }
    } else {
      for (size_t j = min; j < newdeg; j++) {
        res[j] = a[j];
      }
    }
  }
  Polynom<T> result(newdeg, res);
  delete[] res;
  return result;
}

//Вычитание -
template <class T>
Polynom<T> operator-(const Polynom<T>& a, const Polynom<T>& b) {
  size_t newdeg;
  size_t min;
  bool flag = true;
  if (a.Degree() < b.Degree()) {
    newdeg = b.Degree() + 1;
    min = a.Degree() + 1;
  } else {
    newdeg = a.Degree() + 1;
    min = b.Degree() + 1;
    flag = false;
  }
  if (a.Degree() == b.Degree()) {
    size_t i = a.Degree();
    while ((a[i] - b[i] == 0) && (i >= 0) && (i <= a.Degree())) {
      if (i != 0) {
        newdeg--;
        min--;
        i--;
      } else {
        newdeg--;
        min--;
        break;
      }
    }
  }
  T* res;
  if (newdeg == 0) {
    res = nullptr;
  } else {
    res = new T[newdeg];
    for (size_t i = 0; i < min; i++) {
      res[i] = a[i] - b[i];
    }
    if (flag) {
      for (size_t j = min; j < newdeg; j++) {
        res[j] = -b[j];
      }
    } else {
      for (size_t j = min; j < newdeg; j++) {
        res[j] = a[j];
      }
    }
  }
  Polynom<T> result(newdeg, res);
  delete[] res;
  return result;
}

//Умножение на полином *
template <class T>
Polynom<T> operator*(const Polynom<T>& a, const Polynom<T>& b) {
  size_t sa = a.Degree();
  size_t sb = b.Degree();
  size_t newdeg;
  T* res;
  if (((sa == 0) && (a[0] == 0)) || ((sb == 0) && (b[0] == 0))) {
    res = nullptr;
    newdeg = 0;
  } else {
    newdeg = sa + sb + 1;
    res = new T[newdeg];
    for (size_t k = 0; k < newdeg; k++) {
      res[k] = 0;
    }
    for (size_t i = 0; i < sa + 1; i++) {
      for (size_t j = 0; j < sb + 1; j++) {
        res[i + j] += a[i] * b[j];
      }
    }
  }
  Polynom<T> result(newdeg, res);
  delete[] res;
  return result;
}

//Деление на полином /
template <class T>
Polynom<double> operator/(const Polynom<T>& a, const Polynom<T>& b) {
  if (((a.Degree() == 0) && (a[0] == 0)) ||
      ((b.Degree() == 0) && (b[0] == 0))) {
    auto* copyofa = new double[a.Degree() + 1];
    for (size_t i = 0; i < a.Degree() + 1; i++) {
      copyofa[i] = a[i];
    }
    Polynom<double> acp(a.Degree() + 1, copyofa);
    delete[] copyofa;
    return acp;
  } else {
    auto* copyofa = new double[a.Degree() + 1];
    for (size_t i = 0; i < a.Degree() + 1; i++) {
      copyofa[i] = a[i];
    }
    Polynom<double> acp(a.Degree() + 1, copyofa);
    auto* copyofb = new double[b.Degree() + 1];
    for (size_t i = 0; i < b.Degree() + 1; i++) {
      copyofb[i] = b[i];
    }
    Polynom<double> bcp(b.Degree() + 1, copyofb);
    delete[] copyofa;
    delete[] copyofb;
    if (a.Degree() < b.Degree()) {
      double* r;
      r = new double[1];
      r[0] = 0;
      Polynom<double> result(1, r);
      delete[] r;
      return result;
    } else {
      double* res;
      res = new double[a.Degree() - b.Degree() + 1];
      while (acp.Degree() >= bcp.Degree()) {
        size_t delimoe = acp.Degree();
        size_t delitel = bcp.Degree();
        res[delimoe - delitel] = acp[acp.Degree()] / bcp[bcp.Degree()];
        auto* min = new double[delimoe - delitel + 1];
        for (size_t i = 0; i < delimoe - delitel + 1; i++) {
          min[i] = 0;
        }
        min[delimoe - delitel] = res[delimoe - delitel];
        Polynom<double> m(delimoe - delitel + 1, min);
        acp = acp - m * bcp;
        delete[] min;
      }
      Polynom<double> result(a.Degree() - b.Degree() + 1, res);
      delete[] res;
      return result;
    }
  }
}

//Остаток от деления на полином %
Polynom<double> operator%(const Polynom<double>& a, const Polynom<double>& b) {
  if ((b.Degree() == 0) && (b[0] == 0)) {
    auto* r = new double[a.Degree()];
    for (size_t i = 0; i < a.Degree() + 1; i++) {
      r[i] = a[i];
    }
    Polynom<double> res(a.Degree() + 1, r);
    delete[] r;
    return res;
  } else {
    Polynom<double> ch = a / b;
    ch *= b;
    Polynom<double> result = a - ch;
    return result;
  }
}

#endif  // INCLUDE_POLYNOMIAL_HPP_
