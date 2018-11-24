#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <array>
#include <cmath>


template<typename T,size_t N>
struct vec
{
    std::array<T,N> arr;

    T& operator [](size_t i)
    {
        return arr[i];
    }
    const T& operator [](size_t i)const
    {
        return arr[i];
    }
    vec<T,N> operator +(const vec<T,N>& v)const
    {
        vec<T,N> vec;
        for(size_t i = 0; i < N; i++)
            vec[i] = arr[i] + v[i];
        return vec;
    }
    vec<T,N> operator -(const vec<T,N>& v)const
    {
        vec<T,N> vec;
        for(size_t i = 0; i < N; i++)
            vec[i] = arr[i] - v[i];
        return vec;
    }
    vec<T,N> operator *(const T& scalar)const
    {
        vec<T,N> vec;
        for(size_t i = 0; i < N; i++)
            vec[i] = arr[i] * scalar;
        return vec;
    }
    vec<T,N> operator /(const T& scalar)const
    {
        vec<T,N> vec;
        for(size_t i = 0; i < N; i++)
            vec[i] = arr[i] / scalar;
        return vec;
    }
    void operator += (const vec<T,N>& v)
    {
        *this = *this + v;
    }
    void operator -= (const vec<T,N>& v)
    {
        *this = *this - v;
    }
    void operator *= (const T& scalar)
    {
        *this = *this* scalar;
    }
    void operator /= (const T& scalar)
    {
        *this = *this / scalar;
    }
    T dot(const vec<T,N>& v)const
    {
        T sum = 0;
        for(size_t i = 0; i < N; i++)
        {
            sum += arr[i] * v[i];
        }
        return sum;
    }
    T length()const
    {
        return sqrt(this->dot(*this));
    }
    vec<T,N> normalized()const
    {
        return *this / this->length();
    }
    void normalized()
    {
        *this = this->normalized();
    }
};
template<typename T,size_t N>
static std::ostream& operator <<(std::ostream& of, const vec<T,N>& v)
{
    for(size_t i = 0; i < N; i++)
        of<<v[i]<<" ";
    return of;
}
template<typename T>
static vec<T,3> cross(const vec<T,3>& v, const vec<T,3>& w)
{
    vec<T,3> u;
    u[0] = v[1] * w[2] - v[2] * w[1];
    u[1] = v[2] * w[0] - v[0] * w[2];
    u[2] = v[0] * w[1] - v[1] * w[0];
    return u;
}

template<typename T,size_t N>
struct mat
{
    std::array<vec<T,N>,N> arr;

    mat<T,N>& identity()
    {
        memset(arr.data(),0,sizeof(arr));
        for(size_t i = 0; i < N; i++)
            arr[i][i] = 1.0;
        return *this;
    }
    vec<T,N>& operator [](size_t i)
    {
        return arr[i];
    }
    const vec<T,N>& operator [](size_t i)const
    {
        return arr[i];
    }
    mat<T,N> operator +(const mat<T,N>& m)const
    {
        mat<T,N> mat;
        for(size_t i = 0; i < N; i++)
        {
            mat[i] = arr[i] + m[i];
        }
        return mat;
    }
    mat<T,N> operator -(const mat<T,N>& m)const
    {
        mat<T,N> mat;
        for(size_t i = 0; i < N; i++)
        {
            mat[i] = arr[i] - m[i];
        }
        return mat;
    }
    mat<T,N> operator *(const T& scalar)const
    {
        mat<T,N> mat;
        for(size_t i = 0; i < N; i++)
        {
            mat[i] = arr[i] * scalar;
        }
        return mat;
    }
    mat<T,N> operator /(const T& scalar)const
    {
        mat<T,N> mat;
        for(size_t i = 0; i < N; i++)
        {
            mat[i] = arr[i] / scalar;
        }
        return mat;
    } 
    mat<T,N> operator *(const mat<T,N>& m)const
    {
        mat<T,N> matrix;
        mat<T,N> tras = m.trasposta();
        for(size_t i = 0; i < N; i++)
            for(size_t j = 0; j < N; j++)
                matrix[i][j] = arr[i].dot(tras[j]);
        return matrix;
    }
    void operator *=(const mat<T,N>& m)
    {
        *this = *this * m;
    }
    void operator +=(const T& scalar)
    {
        *this = *this + scalar;
    }
    void operator -=(const T& scalar)
    {
        *this = *this - scalar;
    }
    void operator *=(const T& scalar)
    {
        *this = *this * scalar;
    }
    void operator /=(const T& scalar)
    {
        *this = *this / scalar;
    }
    T determinant()const//using gaussian reduction method
    {
        mat<T,N> copy = *this;
        int sign = 1;
        size_t times = N - 1;
        for(size_t i = 0; i < times; i++)
        {
            if(copy[i][i] == 0)
            {
                size_t notNullRow = 0;
                for(size_t j = i + 1; j < N ; j++)
                {
                    if(copy[j][i] != 0)
                    {
                        notNullRow = j;
                        break;
                    }
                }

                if(notNullRow == 0)//all the bottom rows are null
                    return 0;

                std::swap(copy[i],copy[notNullRow]);//swap the rows
                sign *= -1;
            }

            for(size_t j = i + 1; j < N ; j++)
            {
                if(copy[j][i] == 0)
                    continue;
                
                T mul = -(copy[j][i] / copy[i][i]);
                vec<T,N> v = copy[i] * mul;
                copy[j] += v;
            }
            
        }

        T det = sign;
        for(size_t i = 0; i < N; i++)
        {
            det *= copy[i][i];
        }

        return det;
    }
    mat<T,N-1> reduced(size_t pos_i,size_t pos_j)const
    {
        mat<T,N - 1> reduced;
    
        size_t r_i = 0;
        for(size_t i = 0; i < N; i++)
        {
            if(i == pos_i)//skip row
                continue;

            size_t r_j = 0;
            for(size_t j = 0; j < N; j++)
            {
                if(j == pos_j)//skip column
                    continue;
                reduced[r_i][r_j] = arr[i][j];
                r_j++;
            }
            r_i++;
        }

        return reduced;
    }
    T complemento_algebrico(size_t i, size_t j)const
    {
        if(N == 1)
            return 1;
        
        mat<T,N-1> reduced = this->reduced(i,j);
        T det = reduced.determinant();
        if((i + j) % 2 == 1)
            det = - det;
        return det;
    }
    mat<T,N> trasposta()const
    {
        mat<T,N> tras;
        for(size_t i = 0; i < N; i++)
            for(size_t j = 0; j < N; j++)
                tras[i][j] = arr[j][i];
        return tras;
    }
    mat<T,N> aggiunta()const
    {
        mat<T,N> agg;
        for(size_t i = 0; i < N; i++)
            for(size_t j = 0; j < N; j++)
                agg[i][j] = this->complemento_algebrico(i,j);
        return agg;
    }
    mat<T,N> inverse(bool* exist = nullptr)const
    {
        T det = this->determinant();
        if(det == 0)
        {
            if(exist != nullptr)
                *exist = false;
            return {0};
        }
        det = 1.0 / det;
        if(exist != nullptr)
            *exist = true;

        mat<T,N> inverse = this->aggiunta().trasposta();
        inverse = inverse * det;
        return inverse;
    }
    template<typename T1>
    static size_t rango(const mat<T1,1>& m)
    {
        return 1;
    }
    template<typename T1,size_t N1>
    static size_t rango(const mat<T1,N1>& m)
    {
        if(m.determinant() != 0)
            return N1;
        
        size_t maxRango = 0;
        for(size_t i = 0; i < N1; i++) 
            for(size_t j = 0; j < N1; j++)
            {
                size_t r = rango(m.reduced(i,0));
                if(r > maxRango)
                {
                    if(r == N1 - 1)
                        return r;
                    maxRango = r;
                }
            }
    
        return maxRango;
    }
    size_t rango()const
    {
        return rango(*this);
    }
};
template<typename T,size_t N>
static std::ostream& operator <<(std::ostream& of, const mat<T,N>& v)
{
    for(size_t i = 0; i < N; i++)
    {
        of<<v[i];
        if(i != N - 1)
            of<<std::endl;
    }
    return of;
}

#define TYPE double
typedef mat<TYPE,5> mat5;
typedef mat<TYPE,4> mat4;
typedef mat<TYPE,3> mat3;
typedef mat<TYPE,2> mat2;

typedef vec<TYPE,5> vec5;
typedef vec<TYPE,4> vec4;
typedef vec<TYPE,3> vec3;
typedef vec<TYPE,2> vec2;

int main(int argc, char *argv[])
{
 
    mat4 m = { vec4{1,4,0,0},
               vec4{1,8,0,3},
               vec4{1,2,-1,2},
               vec4{5,2,4,1}};
    std::cout<< m.identity()<<std::endl;
    return 0;
}