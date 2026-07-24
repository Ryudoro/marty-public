#pragma once

#include <cassert>
#include <iostream>
#include <exception>
#include <type_traits>

namespace csl {

template <typename T>
class InitSanitizer;

template <typename T>
struct is_init_sanitizer : std::false_type {};

template <typename T>
struct is_init_sanitizer<InitSanitizer<T>> : std::true_type {};

/**
 * @brief Encapsulates a value of a given to ensure that it initialized when
 * used.
 *
 * @details An object of type InitSanitizer<T> can be used as a number for
 * assignment and will be automatically converted to the encapsulated value
 * when necessary. When assigned, the value becomes safe (initially not). When
 * the encapsulated value is taken from the object, it checks that it is indeed
 * initialized before returning it to ensure that there is no garbage result.
 *
 * @tparam T Type of the encapsulated object.
 */
template <typename T>
class InitSanitizer {

  public:
    constexpr InitSanitizer()
    {
    }

    constexpr InitSanitizer(const char * t_name) : name(t_name)
    {
    }

    InitSanitizer(T const &t)
    {
        this->operator=(t);
    }

    InitSanitizer(const char * t_name, T const &t) : name(t_name)
    {
        this->operator=(t);
    }
    
    InitSanitizer &operator=(T const &t)
    {  // Operator already defined before merging "compare"
#ifdef DEBUG
        std::cout << "Calling InitSanitizer &operator=(T const &t) on " << name << "\n";
#endif
        m_safe  = true;
        m_value = t;
        return *this;
    }    

    // ATTENTION: do NOT define the assignment operator in the following comment.
    //            In fact, it will conflict to the one uncommented below 
    //            (i.e. InitSanitizer &operator=(const InitSanitizer<T> &other)).
    //            By defining them both there will be a warning from -Wdeprecated_copy.
    //            Alternatively, if one defines only the one with the two templates, then the other 
    //            is defined by default, and the "name" member WILL be modified.
//     template<typename U>
//     InitSanitizer &operator=(const InitSanitizer<U> &other)
//     { 
// #ifdef DEBUG
//         std::cout << "Calling InitSanitizer &operator=(InitSanitizer<U> const &other) on " << name << "\n";
// #endif
//         this->operator=(other.get());
//         return *this;
//     }
    
    InitSanitizer &operator=(const InitSanitizer<T> &other)
    { 
#ifdef DEBUG
        std::cout << "Calling InitSanitizer &operator=(InitSanitizer<T> const &other) on " << name << "\n";
#endif
        m_safe=other.hasValue();
        m_value= (m_safe ? other.m_value : m_value);
        return *this;
    }
    
    InitSanitizer(InitSanitizer<T> const &other) : 
      name(other.getName()),
      m_value(other.m_value),
      m_safe(other.m_safe)
    {
    }

    bool hasValue() const
    {
        return m_safe;
    }

    T get() const
    {
#ifdef DEBUG
        std::cout << "Called .get()\n";
#endif
        if (!m_safe) {
            std::cerr << "Error: param \"" << getName() << "\" is used ";
            std::cerr << "uninitialized, please assign it a m_value using ";
            std::cerr << "standard m_value assignement.\n";
            throw std::runtime_error("Uninitialized value");
        }
        return m_value;
    }

    void reset()
    {
        m_safe = false;
    }
    
    template <typename U>
    bool operator==(const InitSanitizer<U> &other) const
    {
        if (hasValue() != other.hasValue())
            return false;
        if (!hasValue())
            return true;

        using Common = std::common_type_t<T, U>;
        return static_cast<Common>(get()) == static_cast<Common>(other.get());
    }

    template <typename U>
        requires(!is_init_sanitizer<std::remove_cvref_t<U>>::value
                 && requires(const T &lhs, const U &rhs) { lhs == rhs; })
    bool operator==(const U &other) const
    {
        return hasValue() ? (m_value == other) : false;
    }

    template <typename U>
    bool operator!=(const InitSanitizer<U> &other) const
    {
        return !(*this == other);
    }

    template <typename U>
        requires(!is_init_sanitizer<std::remove_cvref_t<U>>::value
                 && requires(const T &lhs, const U &rhs) { lhs == rhs; })
    bool operator!=(const U &other) const
    {
        return !(*this == other);
    }
    
    operator T() const
    {
        return get();
    }

    std::string getName() const 
    {
      return name;
    }

    void setName(const std::string &newname) 
    {
        name=std::move(newname);
    }
    
    void print(std::ostream &out=std::cout) const
    {
        out << getName() << " = ";
        (hasValue() ?  (out << get()) : (out << "uninitialized")) << '\n'; 
    }
  
  private:
    std::string name="unnamed_var";
    T    m_value;
    bool m_safe{false};
};

} // namespace csl


template <class T>
std::ostream& operator << (std::ostream &out, const csl::InitSanitizer<T> &var)  
{
  var.print(out);
  return out;
}
