// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Mvec.hpp
// This file is part of the Garamon for c3ga2.
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Mvec.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Class to define a multivector and its basic operators in the Geometric algebra of c3ga2.


// Anti-doublon
#ifndef C3GA2_MULTI_VECTOR_HPP__
#define C3GA2_MULTI_VECTOR_HPP__
#pragma once

// External Includes
#include <Eigen/Core>
#include <list>
#include <iostream>
#include <cmath>
#include <limits>

// Internal Includes
#include "c3ga2/Utility.hpp"
#include "c3ga2/Constants.hpp"

#include "c3ga2/Outer.hpp"
#include "c3ga2/Inner.hpp"
#include "c3ga2/Geometric.hpp"

#include "c3ga2/OuterExplicit.hpp"
#include "c3ga2/InnerExplicit.hpp"
#include "c3ga2/GeometricExplicit.hpp"

/*!
 * @namespace c3ga2
 */
namespace c3ga2{


    /// \class Kvec
    /// \brief class defining a single grade component of a multivector.
    template<class T>
    struct Kvec{

        Eigen::Matrix<T, Eigen::Dynamic, 1> vec;  /*!< dynamic vector of Eigen Library */

        unsigned int grade; /*!< grade k of the k-vector */

        /// \brief operator == to test the equality between 2 k-vectors.
        bool operator == (const Kvec<T>& other) const {
            return vec == other.vec;
        }
    };



    /// \class Mvec
    /// \brief class defining multivectors.
    template<typename T = double>
    class Mvec {

    protected:
        std::list<Kvec<T>> mvData;  /*!< set of k-vectors, mapped by grade */
        unsigned int gradeBitmap;   /*!< ith bit to 1 if grade i is contained in the multivector */

    public:

        /// \brief Default constructor, generate an empty multivector equivalent to the scalar 0.
        Mvec();

        /// \brief Copy constructor
        /// \param mv - the multivector to be copied
        Mvec(const Mvec& mv);

        /// \brief Copy constructor of Mvec
        /// \param mv - the multivector which has to be copied
        Mvec(Mvec&& mv); // move constructor


        template <typename U>
        friend class Mvec;

        /// \brief Copy constructor of Mvec from different types
        /// \param mv - the multivector with another template type
        template<typename U>
        Mvec<T>(const Mvec<U> &mv);

        /// \brief Constructor of Mvec from a scalar
        /// \param val - scalar value
        template<typename S>
        Mvec(const S val);

        /// Destructor
        ~Mvec();

        /// \brief Overload the assignment operator. No need any copy when the argument is an R-value, just need to move.
        /// \param mv - Mvec used as a Rvalue
        /// \return assign other to this
        Mvec& operator=(Mvec&& mv);

        /// \brief Overload the assignment operator
        /// \param mv - Mvec
        /// \return assign mv to this object
        Mvec& operator=(const Mvec& mv);

        /// \brief defines the addition between two Mvec
        /// \param mv2 - second operand of type Mvec
        /// \return this + mv2
        Mvec operator+(const Mvec &mv2) const;

        /// \brief defines the addition between a Mvec and a scalar
        /// \param value - second operand (scalar)
        /// \return this + scalar
        template<typename S>
        Mvec operator+(const S &value) const;

        /// \brief defines the addition between a scalar and a Mvec
        /// \param value a scalar
        /// \param mv - second operand of type Mvec
        /// \return scalar + mv
        template<typename U, typename S>
        friend Mvec<U> operator+(const S &value, const Mvec<U> &mv);

        /// \brief Overload the += operator, corresponds to this += mv
        /// \param mv - Mvec to be added to this object
        /// \return this += mv
        Mvec& operator+=(const Mvec& mv);

        /// \brief defines the opposite of a multivector
        /// \param mv: operand of type Mvec
        /// \return -mv
        template<typename U>
        friend Mvec<U> operator-(const Mvec<U> &mv); // unary operator -mv1

        /// \brief defines the difference between two Mvec
        /// \param mv2 - second operand of type Mvec
        /// \return this - mv2
        Mvec<T> operator-(const Mvec<T> &mv2) const;

        /// \brief defines the difference between a Mvec and a scalar
        /// \param value - second operand (scalar)
        /// \return this - scalar
        template<typename S>
        Mvec operator-(const S &value) const;

        /// \brief defines the difference between a scalar and a Mvec
        /// \param value a scalar
        /// \param mv - second operand of type Mvec
        /// \return scalar - mvX
        template<typename U, typename S>
        friend Mvec<U> operator-(const S &value, const Mvec<U> &mv);

        /// \brief Overload the -= operator, corresponds to this -= mv
        /// \param mv - Mvec to be added to this object
        /// \return this -= mv
        Mvec& operator-=(const Mvec& mv);

        /// \brief defines the outer product between two multivectors
        /// \param mv2 - a multivector
        /// \return this^mv2
        Mvec operator^(const Mvec &mv2) const;

        /// \brief defines the outer product a multivector and a scalar
        /// \param value - a scalar
        /// \return this^value
        template<typename S>
        Mvec operator^(const S &value) const;

        /// \brief defines the outer product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return value ^ mv
        template<typename U, typename S>
        friend Mvec<U> operator^(const S &value, const Mvec<U> &mv);

        /// \brief Overload the outer product with operator =, corresponds to this ^= mv
        /// \param mv - Mvec to be wedged to this object
        /// \return this ^= mv
        Mvec& operator^=(const Mvec& mv);

        /// \brief defines the inner product between two multivectors
        /// \param mv2 - a multivector
        /// \return this.mv2
        Mvec operator|(const Mvec &mv2) const;

        /// \brief defines the inner product between a multivector and a scalar
        /// \param value - a scalar
        /// \return this.value = 0 (empty multivector)
        template<typename S>
        Mvec operator|(const S &value) const;

        /// \brief defines the inner product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return mv.value = 0 (empty multivector)
        template<typename U, typename S>
        friend Mvec<U> operator|(const S &value, const Mvec<U> &mv);

        /// \brief Overload the inner product with operator =, corresponds to this |= mv
        /// \param mv - Mvec to be computed in the inner product with this object
        /// \return this |= mv
        Mvec& operator|=(const Mvec& mv);

        /// \brief defines the right contraction between two multivectors
        /// \param mv2 - a multivector
        /// \return the right contraction : $this \\lfloor mv2$
        Mvec operator>(const Mvec &mv2) const;

        /// \brief defines the right contraction between a multivector and a scalar
        /// \param value - a scalar
        /// \return $this \\lfloor value$
        template<typename S>
        Mvec operator>(const S &value) const;

        /// \brief defines the right contraction between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return $value \\lfloor this = 0$ (empty multivector)
        template<typename U, typename S>
        friend Mvec<U> operator>(const S &value, const Mvec<U> &mv);

        /// \brief defines the left contraction between two multivectors
        /// \param mv2 - a multivector
        /// \return the left contraction $this \\rfloor mv2$
        Mvec operator<(const Mvec &mv2) const;

        /// \brief defines the left contraction between a multivector and a scalar
        /// \param value - a scalar
        /// \return $value \\rfloor this$ = 0 (empty multivector)
        template<typename S>
        Mvec operator<(const S &value) const;

        /// \brief defines the left contraction between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return $this \\rfloor value$
        template<typename U, typename S>
        friend Mvec<U> operator<(const S &value, const Mvec<U> &mv);

        /// \brief defines the geometric product between two multivectors
        /// \param mv2 - a multivector
        /// \return mv1*mv2
        Mvec operator*(const Mvec &mv2) const;

    
        /// \brief defines the outer product between two multivectors, where the second multivector is dualized during the product computation.
        /// \param mv2 - a multivector that will be dualized during the wedge with the calling multivector.
        /// \return a multivector.
        Mvec<T> outerPrimalDual(const Mvec<T> &mv2) const;

        /// \brief defines the outer product between two multivectors, where the first multivector is dualized during the product computation.
        /// \param mv2 - a primal form of a multivector; the object multivector will be dualized during the wedge with the calling multivector.
        /// \return a multivector.
        Mvec<T> outerDualPrimal(const Mvec<T> &mv2) const;

        /// \brief defines the outer product between two multivectors, where both multivectors are dualized during the product computation.
        /// \param mv2 - a primal form of a multivector; the object multivector will be dualized during the wedge with the calling multivector.
        /// \return a multivector.
        Mvec<T> outerDualDual(const Mvec<T> &mv2) const;
    
        

        /// \brief defines the scalar product between two multivectors (sum of the inner products between same grade pairs from the 2 multivectors)
        /// \param mv2 - a multivector
        /// \return a scalar
        Mvec<T> scalarProduct(const Mvec<T> &mv2) const;

        /// \brief defines the Hestenes product between two multivectors (Inner product - scalar product)
        /// \param mv2 - a multivector
        /// \return a multivector.
        Mvec<T> hestenesProduct(const Mvec<T> &mv2) const;

        /// \brief defines the geometric product between a multivector and a scalar
        /// \param value - a scalar
        /// \return mv2*value
        template<typename S>
        Mvec operator*(const S &value) const;

        /// \brief defines the geometric product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return value*mv2
        template<typename U, typename S>
        friend Mvec<U> operator*(const S &value, const Mvec<U> &mv);

        /// \brief Overload the geometric product with operator =, corresponds to this *= mv
        /// \param mv - Mvec to be multiplied to this object
        /// \return this *= mv
        Mvec& operator*=(const Mvec& mv);

        /// \brief defines the geometric product with a multivector and the inverse of a second multivector
        /// \param mv2 - a multivector
        /// \return this / mv2
        Mvec operator/(const Mvec &mv2) const;

        /// \brief defines the scalar inverse between a multivector and a scalar
        /// \param value - a scalar
        /// \return this / value
        template<typename S>
        Mvec operator/(const S &value) const;

        /// \brief defines the inverse product between a scalar and a multivector
        /// \param value - a scalar
        /// \param mv - a multivector
        /// \return value / mv
        template<typename U, typename S>
        friend Mvec<U> operator/(const S &value, const Mvec<U> &mv);

        /// \brief Overload the inverse with operator =, corresponds to this /= mv
        /// \param mv - Mvec to be inversed to this object
        /// \return this /= mv
        Mvec& operator/=(const Mvec& mv);

        /// \brief the reverse of a multivector, i.e. if mv = a1^a2^...^an, then reverse(mv) = an^...^a2^a1
        /// \param mv - a multivector
        /// \return reverse of mv
        template<typename U>
        friend Mvec<U> operator~(const Mvec<U> &mv);

        /// \brief the dual of a k-vector is defined as $A_k^* = A_k \\lcont I_n^{-1}$, for a multivector, we just dualize all its components. If the metric is degenerated, this function computes the right complement (mv ^ !mv = I).
        /// \param mv - a multivector
        /// \return the dual of mv
        template<typename U>
        friend Mvec<U> operator!(const Mvec<U> &mv);

        /// \brief boolean operator that tests the equality between two Mvec
        /// \param mv2 - second operand of type Mvec
        /// \return whether two Mvec have the same coefficients
        inline bool operator==(const Mvec& mv2){
            if(gradeBitmap != mv2.gradeBitmap)
                return false;

            return mvData == mv2.mvData;   //// listcaca : ca marche que si les listes sont ordonnees
        }

        /// \brief compute the inverse of a multivector
        /// \return - the inverse of the current multivector
        Mvec<T> inv() const;


            /// \brief operator to test whether two Mvec have not the same coefficients
        /// \param mv2 - second operand of type Mvec
        /// \return boolean that specify the non-equality between two Mvec
        inline bool operator!=(const Mvec& mv2){
            return !(mvData == mv2.mvData);
        }

        /// \brief Display all the non-null basis blades of this objects
        /// \param stream - destination stream
        /// \param mvec - the multivector to be outputed
        /// \return a stream that contains the list of the non-zero element of the multivector
        template<typename U>
        friend std::ostream& operator<< (std::ostream& stream, const Mvec<U> &mvec);

        /// \brief overload the casting operator, using this function, it is now possible to compute : float a = float(mv);
        /// \return the scalar part of the multivector
        operator T () {
            // if no scalar part, return

            if( (gradeBitmap & 1) == 0 )
                return 0;

            // assuming now that the scalar part exists, return it
            return findGrade(0)->vec[0];
        }

        /// \brief Overload the [] operator to assign a basis blade to a multivector. As an example, float a = mv[E12] = 42.
        /// \param idx - the basis vector index related to the query.
        /// \return the coefficient of the multivector corresponding to the "idx" component.
        /// \todo handle the cases when the number of indices is higher than the dimension, and when the index is too high
        T& operator[](const int idx ){

            const unsigned int grade = xorIndexToGrade[idx];
            const unsigned int idxHomogeneous = xorIndexToHomogeneousIndex[idx];

            auto it = mvData.begin();
            while(it != mvData.end()){

                // if grade not reach yet, continue
                if(it->grade < grade)
                    ++it;
                else{

                    // if grade found
                    if(it->grade == grade)
                        return it->vec[idxHomogeneous];

                    // if grade exceed, create it and inster it before the current element
                    Kvec<T> kvec = {Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(binomialArray[grade]),grade};
                    auto it2 = mvData.insert(it,kvec);
                    gradeBitmap |= 1 << (grade);
                    return it2->vec[idxHomogeneous];
                }
            }

            // if the searched element should be added at the end, add it
            Kvec<T> kvec = {Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(binomialArray[grade]),grade};

            auto it2 = mvData.insert(mvData.end(),kvec);
            gradeBitmap |= 1 << (grade);
            return it2->vec[idxHomogeneous];
        }

        /// \brief Overload the [] operator to copy a basis blade of this multivector. As an example, float a = mv[E12].
        /// \param idx - the basis vector index related to the query.
        /// \return the coefficient of the multivector corresponding to the "idx" component.
        const T& operator[](const int idx) const{

            const unsigned int grade = xorIndexToGrade[idx];
            const unsigned int idxHomogeneous = xorIndexToHomogeneousIndex[idx];

            auto it = mvData.begin();
            while(it != mvData.end()){

                // if grade not reach yet, continue
                if(it->grade < grade)
                    ++it;
                else{

                    // if grade found
                    if(it->grade == grade)
                        return it->vec[idxHomogeneous];

                    // if grade exceed, return a reference on zero
                    return zero<T>;
                }
            }

            // searched element not found
            return zero<T>;
        }

/*
        /// \cond DEV
        // Variadic operators
        /// \brief Overload the () operator to assign a basis blade to a multivector. As an example, consider we have a Mvec a and we want to set its e23 component to 4.7, using operator(), it will consists in writing a[2,3]=4.8.
        /// \tparam Arguments - denotes a variadic list of index
        /// \param listIndices - the list of indices
        /// \return the right index in the homogeneous vectorXd
        /// \todo future work + handle the cases when the number of indices is higher than the dimension, and when the index is too high
        template<typename... List>
        T& operator()(List ...listIndices){
            // operation on these arguments
            // First determine the grade, simply use the variadic function sizeof...()
            // This function is correct in terms of indexing
            const int gd = sizeof...(listIndices); //

            // Second, from these parameters, compute the index in the corresponding VectorXd
            const int idx = (Binomial<algebraDimension,gd>::value-1) - computeIdxFromList(algebraDimension,gd,listIndices...);//Binomial<algebraDimension,binomCoef>::value;//binomCoef - sumB(first,listIndices...);//computeOrderingIndex<>(listIndices...);//(Binomial<algebraDimension,gd>::value-1);//-computeOrderingIndex<int,List>(4,2,listIndices...);//computeOrderingIndex<algebraDimension,gd>(listIndices...);//idxVariadic<algebraDimension,gd,first,listIndices...>::value; //computeOrderingIndex<algebraDimension,gd,listIndices...>();//idxVariadic<algebraDimension,gd,listIndices...>::value; // (Binomial<algebraDimension,gd>::value-1)- idxVariadic<algebraDimension,gd,listIndices...>::value

            if(mvData.count(gd)>0)
                return (mvData.at(gd))(idx);

            mvData[gd] = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(Binomial<algebraDimension,gd>::value);

            return mvData[gd](idx);
        }
        /// \endcond
*/

        /// \cond DEV
        /// \brief returns a multivector with one component to 1
        /// \param grade : grade of the component to enable
        /// \param index : index of the parameter in the k-vector (k = grade)
        inline Mvec componentToOne(const unsigned int grade, const int index){
            Kvec<T> kvec;
			kvec.vec=Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(binomialArray[grade]);
			kvec.grade=grade;
            kvec.vec[index] = T(1);
            Mvec mv1;
            mv1.mvData.push_back(kvec);
            mv1.gradeBitmap = 1 << (grade);

            return mv1;
        }
        /// \endcond // do not comment this functions

        /// \cond DEV
        /// \brief create a VectorXd if it has not yet been created
        /// \param grade - grade of the considered kvector
        /// \return nothing
        inline typename std::list<Kvec<T>>::iterator createVectorXdIfDoesNotExist(const unsigned int grade){

            auto it = mvData.begin();
            while(it != mvData.end()){

                // if grade not reach yet, continue
                if(it->grade < grade)
                    ++it;
                else{

                    // if grade found
                    if(it->grade == grade)
                        return it;

                    // if grade exceed, create it and inster it before the current element
                    Kvec<T> kvec;
					kvec.vec=Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(binomialArray[grade]);
					kvec.grade=grade;
                    auto it2 = mvData.insert(it,kvec);
                    gradeBitmap |= 1 << (grade);
                    return it2;
                }
            }

            // if the searched element should be added at the end, add it
            Kvec<T> kvec;
			kvec.vec=Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(binomialArray[grade]);
			kvec.grade=grade;
            auto it2 = mvData.insert(mvData.end(),kvec);
            gradeBitmap |= 1 << (grade);
            return it2;
        }


        /// \cond DEV
        /// \brief modify the element of the multivector whose grade is "grade"
        /// \param grade : the grade to access
        /// \param indexVectorXd : index of the k-vector (with k = "grade")
        /// \return the element of the Mv whose grade is grade and index in the VectorXd is indexVectorXd
        inline T& at(const int grade, const int indexVectorXd){
            auto it = createVectorXdIfDoesNotExist(grade);
            return it->vec.coeffRef(indexVectorXd);  
        }
        /// \endcond // do not comment this functions

        /// \cond DEV
        /// \brief return the element of the multivector whose grade is "grade"
        /// \param grade : the grade to access
        /// \param indexVectorXd : index of the k-vector (with k = "grade")
        /// \return the element of the Mv whose grade is grade and index in the VectorXd is indexVectorXd
        inline T at(const int grade, const int indexVectorXd) const{

            if((gradeBitmap & (1<<(grade))) == 0 )
                return T(0);

            auto it = findGrade(grade);
            return it->vec.coeff(indexVectorXd);
        }
        /// \endcond // do not comment this functions

        /// \brief the L2-norm of the mv is sqrt( abs( mv.mv ) )
        /// \return the L2-norm of the multivector (as a double)
        T inline norm() const {
            return sqrt( fabs( (*this).scalarProduct( this->reverse() ) ));
        };

        /// \brief the L2-norm over 2 of the mv is mv.mv
        /// \return the L2-norm of the multivector (as a double)
        T inline quadraticNorm() const {
            return ( this->reverse() | (*this) );
        };

        /// \brief compute the dual of a multivector (i.e mv* = reverse(mv) * Iinv). If the metric is degenerated, this function computes the right complement (mv ^ !mv = I).
        /// \return - the dual of the multivector
        Mvec<T> dual() const;

        /// \brief compute the reverse of a multivector
        /// \return - the reverse of the multivector
        Mvec<T> reverse() const;

        /// \brief search in the multivector for a Kvec of grade "grade"
        /// \return return a const iterator on the Kvec if exist, else return mvData.end()
        inline typename std::list<Kvec<T>>::const_iterator findGrade(const unsigned int & gradeToFind) const {
            for(auto it = mvData.begin();it != mvData.end();++it){
                if(it->grade==gradeToFind){
                    return it;
                }
            }
            return mvData.end();
        }

        /// \brief search in the multivector for a Kvec of grade "grade"
        /// \return return an iterator on the Kvec if exist, else return mvData.end()
        inline typename std::list<Kvec<T>>::iterator findGrade(const unsigned int & gradeToFind) {
            for(auto it = mvData.begin();it != mvData.end();++it){
                if(it->grade==gradeToFind){
                    return it;
                }
            }
            return mvData.end();
        }

        /// \brief return the (highest) grade of the multivector
        /// \return the highest grade of the multivector
        inline int grade() const {
            return (mvData.begin() == mvData.end()) ? 0 : mvData.rbegin()->grade;
        }

        /// \brief return the all non-zero grades of the multivector (several grades for non-homogeneous multivectors)
        /// \return a list of the grades encountered in the multivector
        std::vector<unsigned int> grades() const;

        /// \brief returns a multivector that contains all the components of this multivector whose grade is i
        /// \return an empty Mvec if the requested element is not part of the multivector, or the multivector that contains only this element if present in the current multivector.
        Mvec grade(const int i) const;

        /// \brief tell whether the multivector has grade component
        /// \param grade - grade of the considered kvector
        /// \return true if multivector has grade component, false else
        inline bool isGrade(const unsigned int grade) const{
            if((grade == 0) && isEmpty()) return true;
            return ( (gradeBitmap & (1<<grade)) != 0 );
        }

        /// \brief partially or completely erase the content of a multivector
        /// \param grade : if < 0, erase all the multivector, else just erase the part of grade "grade".
        void clear(const int grade = -1);

        /// \brief check is a mutivector is empty, i.e. corresponds to 0.
        /// \return True if the multivector is empty, else False.
        inline bool isEmpty() const {
            return mvData.empty();
        }

        /// \brief A multivector is homogeneous if all its components have the same grade.
        /// \return True if the multivector is homogeneous, else False.
        inline bool isHomogeneous() const {
            return mvData.size() < 2; // only one element, or zero element on the list
        }

        /// \brief inplace simplify the multivector such that all the values with a magnitude lower than a epsilon in the Mv are set to 0.
        /// \param epsilon - threshold, with default value the epsilon of the float/double/long double type from numeric_limits.
        void roundZero(const T epsilon = std::numeric_limits<T>::epsilon());

        /// \brief Specify if two multivectors have the same grade.
        /// \param mv - multivector to compare with.
        /// \return true if the two multivectors have the same grade, else return false.
        inline bool sameGrade(const Mvec<T>& mv) const{
            return grade() == mv.grade();
        }

        /// \brief Display the multivector data (per grade value)
        void display() const;

        /// \cond DEV
        /// \brief a function to output all the element of a k-vector in the multivector indexing order.
        /// \tparam T - type of the multivector
        /// \param stream - stream that will contain the result
        /// \param mvec - multivector to be processed
        /// \param gradeMV - the considered grade
        /// \param moreThanOne - true if it the first element to display (should we put a '+' before)
        template<typename U>
        friend void traverseKVector(std::ostream &stream, const Eigen::Matrix<U, Eigen::Dynamic, 1>, unsigned int gradeMV, bool& moreThanOne);
        /// \endcond // do not comment this functions

/*
        /// \cond DEV
        /// \brief functions that enables to  to be able to extract a component of a multivector as a multivector. As an example, mv1.e(1,3) will create a multivector whose 1,3 component will be the component of mv1.
        /// \tparam Arguments - denotes a variadic list of index
        /// \param listIndices - the list of indices
        /// \return a homogeneous multivector filled with zero except the listIndices index
        /// \todo future work + handle the cases when the number of indices is higher than the dimension, and when the index is too high
        template<typename... List>
        Mvec e(List ...listIndices){
            // operation on these arguments
            // First determine the grade, simply use the variadic function sizeof...()
            const int gd = sizeof...(listIndices); //

            // Second, from these parameters, compute the index in the corresponding VectorXd
            const int idx = (Binomial<algebraDimension,gd>::value-1) - computeIdxFromList(algebraDimension,gd,listIndices...);

            Mvec res;
            res.mvData[gd] = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(Binomial<algebraDimension,gd>::value);
            res.mvData[gd](idx) = mvData[gd](idx);

            return res;
        }
        /// \endcond // do not comment this functions
*/


        /// \brief functions that enables to extract one component of this multivector and initialize another mv with this component
        /// \tparam Arguments - denotes a variadic list of index
        /// \param grade : the considered grade
        /// \param sizeOfKVector: C(dimension, grade)
        /// \param indexInKvector: the considerd index
        /// \return a new mv with the right component at the right place
        Mvec extractOneComponent(const int grade, const int sizeOfKVector, const int indexInKvector) const {
            Mvec mv;
            auto itThisMV = this->findGrade(grade);
            if(itThisMV==this->mvData.end()){
                return mv;
            }
            mv = mv.componentToOne(grade,indexInKvector);
            mv.mvData.begin()->vec.coeffRef(indexInKvector) = itThisMV->vec.coeff(indexInKvector);
            return mv;
        }


        /// \brief set of functions that return multivectors containing only the value '1' for the specified component.
        Mvec e01() const {return this->extractOneComponent(1,10, 0);}
        Mvec e1() const {return this->extractOneComponent(1,10, 1);}
        Mvec e2() const {return this->extractOneComponent(1,10, 2);}
        Mvec e3() const {return this->extractOneComponent(1,10, 3);}
        Mvec ei1() const {return this->extractOneComponent(1,10, 4);}
        Mvec e02() const {return this->extractOneComponent(1,10, 5);}
        Mvec e4() const {return this->extractOneComponent(1,10, 6);}
        Mvec e5() const {return this->extractOneComponent(1,10, 7);}
        Mvec e6() const {return this->extractOneComponent(1,10, 8);}
        Mvec ei2() const {return this->extractOneComponent(1,10, 9);}
        Mvec e011() const {return this->extractOneComponent(2,45, 0);}
        Mvec e012() const {return this->extractOneComponent(2,45, 1);}
        Mvec e013() const {return this->extractOneComponent(2,45, 2);}
        Mvec e01i1() const {return this->extractOneComponent(2,45, 3);}
        Mvec e0102() const {return this->extractOneComponent(2,45, 4);}
        Mvec e014() const {return this->extractOneComponent(2,45, 5);}
        Mvec e015() const {return this->extractOneComponent(2,45, 6);}
        Mvec e016() const {return this->extractOneComponent(2,45, 7);}
        Mvec e01i2() const {return this->extractOneComponent(2,45, 8);}
        Mvec e12() const {return this->extractOneComponent(2,45, 9);}
        Mvec e13() const {return this->extractOneComponent(2,45, 10);}
        Mvec e1i1() const {return this->extractOneComponent(2,45, 11);}
        Mvec e102() const {return this->extractOneComponent(2,45, 12);}
        Mvec e14() const {return this->extractOneComponent(2,45, 13);}
        Mvec e15() const {return this->extractOneComponent(2,45, 14);}
        Mvec e16() const {return this->extractOneComponent(2,45, 15);}
        Mvec e1i2() const {return this->extractOneComponent(2,45, 16);}
        Mvec e23() const {return this->extractOneComponent(2,45, 17);}
        Mvec e2i1() const {return this->extractOneComponent(2,45, 18);}
        Mvec e202() const {return this->extractOneComponent(2,45, 19);}
        Mvec e24() const {return this->extractOneComponent(2,45, 20);}
        Mvec e25() const {return this->extractOneComponent(2,45, 21);}
        Mvec e26() const {return this->extractOneComponent(2,45, 22);}
        Mvec e2i2() const {return this->extractOneComponent(2,45, 23);}
        Mvec e3i1() const {return this->extractOneComponent(2,45, 24);}
        Mvec e302() const {return this->extractOneComponent(2,45, 25);}
        Mvec e34() const {return this->extractOneComponent(2,45, 26);}
        Mvec e35() const {return this->extractOneComponent(2,45, 27);}
        Mvec e36() const {return this->extractOneComponent(2,45, 28);}
        Mvec e3i2() const {return this->extractOneComponent(2,45, 29);}
        Mvec ei102() const {return this->extractOneComponent(2,45, 30);}
        Mvec ei14() const {return this->extractOneComponent(2,45, 31);}
        Mvec ei15() const {return this->extractOneComponent(2,45, 32);}
        Mvec ei16() const {return this->extractOneComponent(2,45, 33);}
        Mvec ei1i2() const {return this->extractOneComponent(2,45, 34);}
        Mvec e024() const {return this->extractOneComponent(2,45, 35);}
        Mvec e025() const {return this->extractOneComponent(2,45, 36);}
        Mvec e026() const {return this->extractOneComponent(2,45, 37);}
        Mvec e02i2() const {return this->extractOneComponent(2,45, 38);}
        Mvec e45() const {return this->extractOneComponent(2,45, 39);}
        Mvec e46() const {return this->extractOneComponent(2,45, 40);}
        Mvec e4i2() const {return this->extractOneComponent(2,45, 41);}
        Mvec e56() const {return this->extractOneComponent(2,45, 42);}
        Mvec e5i2() const {return this->extractOneComponent(2,45, 43);}
        Mvec e6i2() const {return this->extractOneComponent(2,45, 44);}
        Mvec e0112() const {return this->extractOneComponent(3,120, 0);}
        Mvec e0113() const {return this->extractOneComponent(3,120, 1);}
        Mvec e011i1() const {return this->extractOneComponent(3,120, 2);}
        Mvec e01102() const {return this->extractOneComponent(3,120, 3);}
        Mvec e0114() const {return this->extractOneComponent(3,120, 4);}
        Mvec e0115() const {return this->extractOneComponent(3,120, 5);}
        Mvec e0116() const {return this->extractOneComponent(3,120, 6);}
        Mvec e011i2() const {return this->extractOneComponent(3,120, 7);}
        Mvec e0123() const {return this->extractOneComponent(3,120, 8);}
        Mvec e012i1() const {return this->extractOneComponent(3,120, 9);}
        Mvec e01202() const {return this->extractOneComponent(3,120, 10);}
        Mvec e0124() const {return this->extractOneComponent(3,120, 11);}
        Mvec e0125() const {return this->extractOneComponent(3,120, 12);}
        Mvec e0126() const {return this->extractOneComponent(3,120, 13);}
        Mvec e012i2() const {return this->extractOneComponent(3,120, 14);}
        Mvec e013i1() const {return this->extractOneComponent(3,120, 15);}
        Mvec e01302() const {return this->extractOneComponent(3,120, 16);}
        Mvec e0134() const {return this->extractOneComponent(3,120, 17);}
        Mvec e0135() const {return this->extractOneComponent(3,120, 18);}
        Mvec e0136() const {return this->extractOneComponent(3,120, 19);}
        Mvec e013i2() const {return this->extractOneComponent(3,120, 20);}
        Mvec e01i102() const {return this->extractOneComponent(3,120, 21);}
        Mvec e01i14() const {return this->extractOneComponent(3,120, 22);}
        Mvec e01i15() const {return this->extractOneComponent(3,120, 23);}
        Mvec e01i16() const {return this->extractOneComponent(3,120, 24);}
        Mvec e01i1i2() const {return this->extractOneComponent(3,120, 25);}
        Mvec e01024() const {return this->extractOneComponent(3,120, 26);}
        Mvec e01025() const {return this->extractOneComponent(3,120, 27);}
        Mvec e01026() const {return this->extractOneComponent(3,120, 28);}
        Mvec e0102i2() const {return this->extractOneComponent(3,120, 29);}
        Mvec e0145() const {return this->extractOneComponent(3,120, 30);}
        Mvec e0146() const {return this->extractOneComponent(3,120, 31);}
        Mvec e014i2() const {return this->extractOneComponent(3,120, 32);}
        Mvec e0156() const {return this->extractOneComponent(3,120, 33);}
        Mvec e015i2() const {return this->extractOneComponent(3,120, 34);}
        Mvec e016i2() const {return this->extractOneComponent(3,120, 35);}
        Mvec e123() const {return this->extractOneComponent(3,120, 36);}
        Mvec e12i1() const {return this->extractOneComponent(3,120, 37);}
        Mvec e1202() const {return this->extractOneComponent(3,120, 38);}
        Mvec e124() const {return this->extractOneComponent(3,120, 39);}
        Mvec e125() const {return this->extractOneComponent(3,120, 40);}
        Mvec e126() const {return this->extractOneComponent(3,120, 41);}
        Mvec e12i2() const {return this->extractOneComponent(3,120, 42);}
        Mvec e13i1() const {return this->extractOneComponent(3,120, 43);}
        Mvec e1302() const {return this->extractOneComponent(3,120, 44);}
        Mvec e134() const {return this->extractOneComponent(3,120, 45);}
        Mvec e135() const {return this->extractOneComponent(3,120, 46);}
        Mvec e136() const {return this->extractOneComponent(3,120, 47);}
        Mvec e13i2() const {return this->extractOneComponent(3,120, 48);}
        Mvec e1i102() const {return this->extractOneComponent(3,120, 49);}
        Mvec e1i14() const {return this->extractOneComponent(3,120, 50);}
        Mvec e1i15() const {return this->extractOneComponent(3,120, 51);}
        Mvec e1i16() const {return this->extractOneComponent(3,120, 52);}
        Mvec e1i1i2() const {return this->extractOneComponent(3,120, 53);}
        Mvec e1024() const {return this->extractOneComponent(3,120, 54);}
        Mvec e1025() const {return this->extractOneComponent(3,120, 55);}
        Mvec e1026() const {return this->extractOneComponent(3,120, 56);}
        Mvec e102i2() const {return this->extractOneComponent(3,120, 57);}
        Mvec e145() const {return this->extractOneComponent(3,120, 58);}
        Mvec e146() const {return this->extractOneComponent(3,120, 59);}
        Mvec e14i2() const {return this->extractOneComponent(3,120, 60);}
        Mvec e156() const {return this->extractOneComponent(3,120, 61);}
        Mvec e15i2() const {return this->extractOneComponent(3,120, 62);}
        Mvec e16i2() const {return this->extractOneComponent(3,120, 63);}
        Mvec e23i1() const {return this->extractOneComponent(3,120, 64);}
        Mvec e2302() const {return this->extractOneComponent(3,120, 65);}
        Mvec e234() const {return this->extractOneComponent(3,120, 66);}
        Mvec e235() const {return this->extractOneComponent(3,120, 67);}
        Mvec e236() const {return this->extractOneComponent(3,120, 68);}
        Mvec e23i2() const {return this->extractOneComponent(3,120, 69);}
        Mvec e2i102() const {return this->extractOneComponent(3,120, 70);}
        Mvec e2i14() const {return this->extractOneComponent(3,120, 71);}
        Mvec e2i15() const {return this->extractOneComponent(3,120, 72);}
        Mvec e2i16() const {return this->extractOneComponent(3,120, 73);}
        Mvec e2i1i2() const {return this->extractOneComponent(3,120, 74);}
        Mvec e2024() const {return this->extractOneComponent(3,120, 75);}
        Mvec e2025() const {return this->extractOneComponent(3,120, 76);}
        Mvec e2026() const {return this->extractOneComponent(3,120, 77);}
        Mvec e202i2() const {return this->extractOneComponent(3,120, 78);}
        Mvec e245() const {return this->extractOneComponent(3,120, 79);}
        Mvec e246() const {return this->extractOneComponent(3,120, 80);}
        Mvec e24i2() const {return this->extractOneComponent(3,120, 81);}
        Mvec e256() const {return this->extractOneComponent(3,120, 82);}
        Mvec e25i2() const {return this->extractOneComponent(3,120, 83);}
        Mvec e26i2() const {return this->extractOneComponent(3,120, 84);}
        Mvec e3i102() const {return this->extractOneComponent(3,120, 85);}
        Mvec e3i14() const {return this->extractOneComponent(3,120, 86);}
        Mvec e3i15() const {return this->extractOneComponent(3,120, 87);}
        Mvec e3i16() const {return this->extractOneComponent(3,120, 88);}
        Mvec e3i1i2() const {return this->extractOneComponent(3,120, 89);}
        Mvec e3024() const {return this->extractOneComponent(3,120, 90);}
        Mvec e3025() const {return this->extractOneComponent(3,120, 91);}
        Mvec e3026() const {return this->extractOneComponent(3,120, 92);}
        Mvec e302i2() const {return this->extractOneComponent(3,120, 93);}
        Mvec e345() const {return this->extractOneComponent(3,120, 94);}
        Mvec e346() const {return this->extractOneComponent(3,120, 95);}
        Mvec e34i2() const {return this->extractOneComponent(3,120, 96);}
        Mvec e356() const {return this->extractOneComponent(3,120, 97);}
        Mvec e35i2() const {return this->extractOneComponent(3,120, 98);}
        Mvec e36i2() const {return this->extractOneComponent(3,120, 99);}
        Mvec ei1024() const {return this->extractOneComponent(3,120, 100);}
        Mvec ei1025() const {return this->extractOneComponent(3,120, 101);}
        Mvec ei1026() const {return this->extractOneComponent(3,120, 102);}
        Mvec ei102i2() const {return this->extractOneComponent(3,120, 103);}
        Mvec ei145() const {return this->extractOneComponent(3,120, 104);}
        Mvec ei146() const {return this->extractOneComponent(3,120, 105);}
        Mvec ei14i2() const {return this->extractOneComponent(3,120, 106);}
        Mvec ei156() const {return this->extractOneComponent(3,120, 107);}
        Mvec ei15i2() const {return this->extractOneComponent(3,120, 108);}
        Mvec ei16i2() const {return this->extractOneComponent(3,120, 109);}
        Mvec e0245() const {return this->extractOneComponent(3,120, 110);}
        Mvec e0246() const {return this->extractOneComponent(3,120, 111);}
        Mvec e024i2() const {return this->extractOneComponent(3,120, 112);}
        Mvec e0256() const {return this->extractOneComponent(3,120, 113);}
        Mvec e025i2() const {return this->extractOneComponent(3,120, 114);}
        Mvec e026i2() const {return this->extractOneComponent(3,120, 115);}
        Mvec e456() const {return this->extractOneComponent(3,120, 116);}
        Mvec e45i2() const {return this->extractOneComponent(3,120, 117);}
        Mvec e46i2() const {return this->extractOneComponent(3,120, 118);}
        Mvec e56i2() const {return this->extractOneComponent(3,120, 119);}
        Mvec e01123() const {return this->extractOneComponent(4,210, 0);}
        Mvec e0112i1() const {return this->extractOneComponent(4,210, 1);}
        Mvec e011202() const {return this->extractOneComponent(4,210, 2);}
        Mvec e01124() const {return this->extractOneComponent(4,210, 3);}
        Mvec e01125() const {return this->extractOneComponent(4,210, 4);}
        Mvec e01126() const {return this->extractOneComponent(4,210, 5);}
        Mvec e0112i2() const {return this->extractOneComponent(4,210, 6);}
        Mvec e0113i1() const {return this->extractOneComponent(4,210, 7);}
        Mvec e011302() const {return this->extractOneComponent(4,210, 8);}
        Mvec e01134() const {return this->extractOneComponent(4,210, 9);}
        Mvec e01135() const {return this->extractOneComponent(4,210, 10);}
        Mvec e01136() const {return this->extractOneComponent(4,210, 11);}
        Mvec e0113i2() const {return this->extractOneComponent(4,210, 12);}
        Mvec e011i102() const {return this->extractOneComponent(4,210, 13);}
        Mvec e011i14() const {return this->extractOneComponent(4,210, 14);}
        Mvec e011i15() const {return this->extractOneComponent(4,210, 15);}
        Mvec e011i16() const {return this->extractOneComponent(4,210, 16);}
        Mvec e011i1i2() const {return this->extractOneComponent(4,210, 17);}
        Mvec e011024() const {return this->extractOneComponent(4,210, 18);}
        Mvec e011025() const {return this->extractOneComponent(4,210, 19);}
        Mvec e011026() const {return this->extractOneComponent(4,210, 20);}
        Mvec e01102i2() const {return this->extractOneComponent(4,210, 21);}
        Mvec e01145() const {return this->extractOneComponent(4,210, 22);}
        Mvec e01146() const {return this->extractOneComponent(4,210, 23);}
        Mvec e0114i2() const {return this->extractOneComponent(4,210, 24);}
        Mvec e01156() const {return this->extractOneComponent(4,210, 25);}
        Mvec e0115i2() const {return this->extractOneComponent(4,210, 26);}
        Mvec e0116i2() const {return this->extractOneComponent(4,210, 27);}
        Mvec e0123i1() const {return this->extractOneComponent(4,210, 28);}
        Mvec e012302() const {return this->extractOneComponent(4,210, 29);}
        Mvec e01234() const {return this->extractOneComponent(4,210, 30);}
        Mvec e01235() const {return this->extractOneComponent(4,210, 31);}
        Mvec e01236() const {return this->extractOneComponent(4,210, 32);}
        Mvec e0123i2() const {return this->extractOneComponent(4,210, 33);}
        Mvec e012i102() const {return this->extractOneComponent(4,210, 34);}
        Mvec e012i14() const {return this->extractOneComponent(4,210, 35);}
        Mvec e012i15() const {return this->extractOneComponent(4,210, 36);}
        Mvec e012i16() const {return this->extractOneComponent(4,210, 37);}
        Mvec e012i1i2() const {return this->extractOneComponent(4,210, 38);}
        Mvec e012024() const {return this->extractOneComponent(4,210, 39);}
        Mvec e012025() const {return this->extractOneComponent(4,210, 40);}
        Mvec e012026() const {return this->extractOneComponent(4,210, 41);}
        Mvec e01202i2() const {return this->extractOneComponent(4,210, 42);}
        Mvec e01245() const {return this->extractOneComponent(4,210, 43);}
        Mvec e01246() const {return this->extractOneComponent(4,210, 44);}
        Mvec e0124i2() const {return this->extractOneComponent(4,210, 45);}
        Mvec e01256() const {return this->extractOneComponent(4,210, 46);}
        Mvec e0125i2() const {return this->extractOneComponent(4,210, 47);}
        Mvec e0126i2() const {return this->extractOneComponent(4,210, 48);}
        Mvec e013i102() const {return this->extractOneComponent(4,210, 49);}
        Mvec e013i14() const {return this->extractOneComponent(4,210, 50);}
        Mvec e013i15() const {return this->extractOneComponent(4,210, 51);}
        Mvec e013i16() const {return this->extractOneComponent(4,210, 52);}
        Mvec e013i1i2() const {return this->extractOneComponent(4,210, 53);}
        Mvec e013024() const {return this->extractOneComponent(4,210, 54);}
        Mvec e013025() const {return this->extractOneComponent(4,210, 55);}
        Mvec e013026() const {return this->extractOneComponent(4,210, 56);}
        Mvec e01302i2() const {return this->extractOneComponent(4,210, 57);}
        Mvec e01345() const {return this->extractOneComponent(4,210, 58);}
        Mvec e01346() const {return this->extractOneComponent(4,210, 59);}
        Mvec e0134i2() const {return this->extractOneComponent(4,210, 60);}
        Mvec e01356() const {return this->extractOneComponent(4,210, 61);}
        Mvec e0135i2() const {return this->extractOneComponent(4,210, 62);}
        Mvec e0136i2() const {return this->extractOneComponent(4,210, 63);}
        Mvec e01i1024() const {return this->extractOneComponent(4,210, 64);}
        Mvec e01i1025() const {return this->extractOneComponent(4,210, 65);}
        Mvec e01i1026() const {return this->extractOneComponent(4,210, 66);}
        Mvec e01i102i2() const {return this->extractOneComponent(4,210, 67);}
        Mvec e01i145() const {return this->extractOneComponent(4,210, 68);}
        Mvec e01i146() const {return this->extractOneComponent(4,210, 69);}
        Mvec e01i14i2() const {return this->extractOneComponent(4,210, 70);}
        Mvec e01i156() const {return this->extractOneComponent(4,210, 71);}
        Mvec e01i15i2() const {return this->extractOneComponent(4,210, 72);}
        Mvec e01i16i2() const {return this->extractOneComponent(4,210, 73);}
        Mvec e010245() const {return this->extractOneComponent(4,210, 74);}
        Mvec e010246() const {return this->extractOneComponent(4,210, 75);}
        Mvec e01024i2() const {return this->extractOneComponent(4,210, 76);}
        Mvec e010256() const {return this->extractOneComponent(4,210, 77);}
        Mvec e01025i2() const {return this->extractOneComponent(4,210, 78);}
        Mvec e01026i2() const {return this->extractOneComponent(4,210, 79);}
        Mvec e01456() const {return this->extractOneComponent(4,210, 80);}
        Mvec e0145i2() const {return this->extractOneComponent(4,210, 81);}
        Mvec e0146i2() const {return this->extractOneComponent(4,210, 82);}
        Mvec e0156i2() const {return this->extractOneComponent(4,210, 83);}
        Mvec e123i1() const {return this->extractOneComponent(4,210, 84);}
        Mvec e12302() const {return this->extractOneComponent(4,210, 85);}
        Mvec e1234() const {return this->extractOneComponent(4,210, 86);}
        Mvec e1235() const {return this->extractOneComponent(4,210, 87);}
        Mvec e1236() const {return this->extractOneComponent(4,210, 88);}
        Mvec e123i2() const {return this->extractOneComponent(4,210, 89);}
        Mvec e12i102() const {return this->extractOneComponent(4,210, 90);}
        Mvec e12i14() const {return this->extractOneComponent(4,210, 91);}
        Mvec e12i15() const {return this->extractOneComponent(4,210, 92);}
        Mvec e12i16() const {return this->extractOneComponent(4,210, 93);}
        Mvec e12i1i2() const {return this->extractOneComponent(4,210, 94);}
        Mvec e12024() const {return this->extractOneComponent(4,210, 95);}
        Mvec e12025() const {return this->extractOneComponent(4,210, 96);}
        Mvec e12026() const {return this->extractOneComponent(4,210, 97);}
        Mvec e1202i2() const {return this->extractOneComponent(4,210, 98);}
        Mvec e1245() const {return this->extractOneComponent(4,210, 99);}
        Mvec e1246() const {return this->extractOneComponent(4,210, 100);}
        Mvec e124i2() const {return this->extractOneComponent(4,210, 101);}
        Mvec e1256() const {return this->extractOneComponent(4,210, 102);}
        Mvec e125i2() const {return this->extractOneComponent(4,210, 103);}
        Mvec e126i2() const {return this->extractOneComponent(4,210, 104);}
        Mvec e13i102() const {return this->extractOneComponent(4,210, 105);}
        Mvec e13i14() const {return this->extractOneComponent(4,210, 106);}
        Mvec e13i15() const {return this->extractOneComponent(4,210, 107);}
        Mvec e13i16() const {return this->extractOneComponent(4,210, 108);}
        Mvec e13i1i2() const {return this->extractOneComponent(4,210, 109);}
        Mvec e13024() const {return this->extractOneComponent(4,210, 110);}
        Mvec e13025() const {return this->extractOneComponent(4,210, 111);}
        Mvec e13026() const {return this->extractOneComponent(4,210, 112);}
        Mvec e1302i2() const {return this->extractOneComponent(4,210, 113);}
        Mvec e1345() const {return this->extractOneComponent(4,210, 114);}
        Mvec e1346() const {return this->extractOneComponent(4,210, 115);}
        Mvec e134i2() const {return this->extractOneComponent(4,210, 116);}
        Mvec e1356() const {return this->extractOneComponent(4,210, 117);}
        Mvec e135i2() const {return this->extractOneComponent(4,210, 118);}
        Mvec e136i2() const {return this->extractOneComponent(4,210, 119);}
        Mvec e1i1024() const {return this->extractOneComponent(4,210, 120);}
        Mvec e1i1025() const {return this->extractOneComponent(4,210, 121);}
        Mvec e1i1026() const {return this->extractOneComponent(4,210, 122);}
        Mvec e1i102i2() const {return this->extractOneComponent(4,210, 123);}
        Mvec e1i145() const {return this->extractOneComponent(4,210, 124);}
        Mvec e1i146() const {return this->extractOneComponent(4,210, 125);}
        Mvec e1i14i2() const {return this->extractOneComponent(4,210, 126);}
        Mvec e1i156() const {return this->extractOneComponent(4,210, 127);}
        Mvec e1i15i2() const {return this->extractOneComponent(4,210, 128);}
        Mvec e1i16i2() const {return this->extractOneComponent(4,210, 129);}
        Mvec e10245() const {return this->extractOneComponent(4,210, 130);}
        Mvec e10246() const {return this->extractOneComponent(4,210, 131);}
        Mvec e1024i2() const {return this->extractOneComponent(4,210, 132);}
        Mvec e10256() const {return this->extractOneComponent(4,210, 133);}
        Mvec e1025i2() const {return this->extractOneComponent(4,210, 134);}
        Mvec e1026i2() const {return this->extractOneComponent(4,210, 135);}
        Mvec e1456() const {return this->extractOneComponent(4,210, 136);}
        Mvec e145i2() const {return this->extractOneComponent(4,210, 137);}
        Mvec e146i2() const {return this->extractOneComponent(4,210, 138);}
        Mvec e156i2() const {return this->extractOneComponent(4,210, 139);}
        Mvec e23i102() const {return this->extractOneComponent(4,210, 140);}
        Mvec e23i14() const {return this->extractOneComponent(4,210, 141);}
        Mvec e23i15() const {return this->extractOneComponent(4,210, 142);}
        Mvec e23i16() const {return this->extractOneComponent(4,210, 143);}
        Mvec e23i1i2() const {return this->extractOneComponent(4,210, 144);}
        Mvec e23024() const {return this->extractOneComponent(4,210, 145);}
        Mvec e23025() const {return this->extractOneComponent(4,210, 146);}
        Mvec e23026() const {return this->extractOneComponent(4,210, 147);}
        Mvec e2302i2() const {return this->extractOneComponent(4,210, 148);}
        Mvec e2345() const {return this->extractOneComponent(4,210, 149);}
        Mvec e2346() const {return this->extractOneComponent(4,210, 150);}
        Mvec e234i2() const {return this->extractOneComponent(4,210, 151);}
        Mvec e2356() const {return this->extractOneComponent(4,210, 152);}
        Mvec e235i2() const {return this->extractOneComponent(4,210, 153);}
        Mvec e236i2() const {return this->extractOneComponent(4,210, 154);}
        Mvec e2i1024() const {return this->extractOneComponent(4,210, 155);}
        Mvec e2i1025() const {return this->extractOneComponent(4,210, 156);}
        Mvec e2i1026() const {return this->extractOneComponent(4,210, 157);}
        Mvec e2i102i2() const {return this->extractOneComponent(4,210, 158);}
        Mvec e2i145() const {return this->extractOneComponent(4,210, 159);}
        Mvec e2i146() const {return this->extractOneComponent(4,210, 160);}
        Mvec e2i14i2() const {return this->extractOneComponent(4,210, 161);}
        Mvec e2i156() const {return this->extractOneComponent(4,210, 162);}
        Mvec e2i15i2() const {return this->extractOneComponent(4,210, 163);}
        Mvec e2i16i2() const {return this->extractOneComponent(4,210, 164);}
        Mvec e20245() const {return this->extractOneComponent(4,210, 165);}
        Mvec e20246() const {return this->extractOneComponent(4,210, 166);}
        Mvec e2024i2() const {return this->extractOneComponent(4,210, 167);}
        Mvec e20256() const {return this->extractOneComponent(4,210, 168);}
        Mvec e2025i2() const {return this->extractOneComponent(4,210, 169);}
        Mvec e2026i2() const {return this->extractOneComponent(4,210, 170);}
        Mvec e2456() const {return this->extractOneComponent(4,210, 171);}
        Mvec e245i2() const {return this->extractOneComponent(4,210, 172);}
        Mvec e246i2() const {return this->extractOneComponent(4,210, 173);}
        Mvec e256i2() const {return this->extractOneComponent(4,210, 174);}
        Mvec e3i1024() const {return this->extractOneComponent(4,210, 175);}
        Mvec e3i1025() const {return this->extractOneComponent(4,210, 176);}
        Mvec e3i1026() const {return this->extractOneComponent(4,210, 177);}
        Mvec e3i102i2() const {return this->extractOneComponent(4,210, 178);}
        Mvec e3i145() const {return this->extractOneComponent(4,210, 179);}
        Mvec e3i146() const {return this->extractOneComponent(4,210, 180);}
        Mvec e3i14i2() const {return this->extractOneComponent(4,210, 181);}
        Mvec e3i156() const {return this->extractOneComponent(4,210, 182);}
        Mvec e3i15i2() const {return this->extractOneComponent(4,210, 183);}
        Mvec e3i16i2() const {return this->extractOneComponent(4,210, 184);}
        Mvec e30245() const {return this->extractOneComponent(4,210, 185);}
        Mvec e30246() const {return this->extractOneComponent(4,210, 186);}
        Mvec e3024i2() const {return this->extractOneComponent(4,210, 187);}
        Mvec e30256() const {return this->extractOneComponent(4,210, 188);}
        Mvec e3025i2() const {return this->extractOneComponent(4,210, 189);}
        Mvec e3026i2() const {return this->extractOneComponent(4,210, 190);}
        Mvec e3456() const {return this->extractOneComponent(4,210, 191);}
        Mvec e345i2() const {return this->extractOneComponent(4,210, 192);}
        Mvec e346i2() const {return this->extractOneComponent(4,210, 193);}
        Mvec e356i2() const {return this->extractOneComponent(4,210, 194);}
        Mvec ei10245() const {return this->extractOneComponent(4,210, 195);}
        Mvec ei10246() const {return this->extractOneComponent(4,210, 196);}
        Mvec ei1024i2() const {return this->extractOneComponent(4,210, 197);}
        Mvec ei10256() const {return this->extractOneComponent(4,210, 198);}
        Mvec ei1025i2() const {return this->extractOneComponent(4,210, 199);}
        Mvec ei1026i2() const {return this->extractOneComponent(4,210, 200);}
        Mvec ei1456() const {return this->extractOneComponent(4,210, 201);}
        Mvec ei145i2() const {return this->extractOneComponent(4,210, 202);}
        Mvec ei146i2() const {return this->extractOneComponent(4,210, 203);}
        Mvec ei156i2() const {return this->extractOneComponent(4,210, 204);}
        Mvec e02456() const {return this->extractOneComponent(4,210, 205);}
        Mvec e0245i2() const {return this->extractOneComponent(4,210, 206);}
        Mvec e0246i2() const {return this->extractOneComponent(4,210, 207);}
        Mvec e0256i2() const {return this->extractOneComponent(4,210, 208);}
        Mvec e456i2() const {return this->extractOneComponent(4,210, 209);}
        Mvec e01123i1() const {return this->extractOneComponent(5,252, 0);}
        Mvec e0112302() const {return this->extractOneComponent(5,252, 1);}
        Mvec e011234() const {return this->extractOneComponent(5,252, 2);}
        Mvec e011235() const {return this->extractOneComponent(5,252, 3);}
        Mvec e011236() const {return this->extractOneComponent(5,252, 4);}
        Mvec e01123i2() const {return this->extractOneComponent(5,252, 5);}
        Mvec e0112i102() const {return this->extractOneComponent(5,252, 6);}
        Mvec e0112i14() const {return this->extractOneComponent(5,252, 7);}
        Mvec e0112i15() const {return this->extractOneComponent(5,252, 8);}
        Mvec e0112i16() const {return this->extractOneComponent(5,252, 9);}
        Mvec e0112i1i2() const {return this->extractOneComponent(5,252, 10);}
        Mvec e0112024() const {return this->extractOneComponent(5,252, 11);}
        Mvec e0112025() const {return this->extractOneComponent(5,252, 12);}
        Mvec e0112026() const {return this->extractOneComponent(5,252, 13);}
        Mvec e011202i2() const {return this->extractOneComponent(5,252, 14);}
        Mvec e011245() const {return this->extractOneComponent(5,252, 15);}
        Mvec e011246() const {return this->extractOneComponent(5,252, 16);}
        Mvec e01124i2() const {return this->extractOneComponent(5,252, 17);}
        Mvec e011256() const {return this->extractOneComponent(5,252, 18);}
        Mvec e01125i2() const {return this->extractOneComponent(5,252, 19);}
        Mvec e01126i2() const {return this->extractOneComponent(5,252, 20);}
        Mvec e0113i102() const {return this->extractOneComponent(5,252, 21);}
        Mvec e0113i14() const {return this->extractOneComponent(5,252, 22);}
        Mvec e0113i15() const {return this->extractOneComponent(5,252, 23);}
        Mvec e0113i16() const {return this->extractOneComponent(5,252, 24);}
        Mvec e0113i1i2() const {return this->extractOneComponent(5,252, 25);}
        Mvec e0113024() const {return this->extractOneComponent(5,252, 26);}
        Mvec e0113025() const {return this->extractOneComponent(5,252, 27);}
        Mvec e0113026() const {return this->extractOneComponent(5,252, 28);}
        Mvec e011302i2() const {return this->extractOneComponent(5,252, 29);}
        Mvec e011345() const {return this->extractOneComponent(5,252, 30);}
        Mvec e011346() const {return this->extractOneComponent(5,252, 31);}
        Mvec e01134i2() const {return this->extractOneComponent(5,252, 32);}
        Mvec e011356() const {return this->extractOneComponent(5,252, 33);}
        Mvec e01135i2() const {return this->extractOneComponent(5,252, 34);}
        Mvec e01136i2() const {return this->extractOneComponent(5,252, 35);}
        Mvec e011i1024() const {return this->extractOneComponent(5,252, 36);}
        Mvec e011i1025() const {return this->extractOneComponent(5,252, 37);}
        Mvec e011i1026() const {return this->extractOneComponent(5,252, 38);}
        Mvec e011i102i2() const {return this->extractOneComponent(5,252, 39);}
        Mvec e011i145() const {return this->extractOneComponent(5,252, 40);}
        Mvec e011i146() const {return this->extractOneComponent(5,252, 41);}
        Mvec e011i14i2() const {return this->extractOneComponent(5,252, 42);}
        Mvec e011i156() const {return this->extractOneComponent(5,252, 43);}
        Mvec e011i15i2() const {return this->extractOneComponent(5,252, 44);}
        Mvec e011i16i2() const {return this->extractOneComponent(5,252, 45);}
        Mvec e0110245() const {return this->extractOneComponent(5,252, 46);}
        Mvec e0110246() const {return this->extractOneComponent(5,252, 47);}
        Mvec e011024i2() const {return this->extractOneComponent(5,252, 48);}
        Mvec e0110256() const {return this->extractOneComponent(5,252, 49);}
        Mvec e011025i2() const {return this->extractOneComponent(5,252, 50);}
        Mvec e011026i2() const {return this->extractOneComponent(5,252, 51);}
        Mvec e011456() const {return this->extractOneComponent(5,252, 52);}
        Mvec e01145i2() const {return this->extractOneComponent(5,252, 53);}
        Mvec e01146i2() const {return this->extractOneComponent(5,252, 54);}
        Mvec e01156i2() const {return this->extractOneComponent(5,252, 55);}
        Mvec e0123i102() const {return this->extractOneComponent(5,252, 56);}
        Mvec e0123i14() const {return this->extractOneComponent(5,252, 57);}
        Mvec e0123i15() const {return this->extractOneComponent(5,252, 58);}
        Mvec e0123i16() const {return this->extractOneComponent(5,252, 59);}
        Mvec e0123i1i2() const {return this->extractOneComponent(5,252, 60);}
        Mvec e0123024() const {return this->extractOneComponent(5,252, 61);}
        Mvec e0123025() const {return this->extractOneComponent(5,252, 62);}
        Mvec e0123026() const {return this->extractOneComponent(5,252, 63);}
        Mvec e012302i2() const {return this->extractOneComponent(5,252, 64);}
        Mvec e012345() const {return this->extractOneComponent(5,252, 65);}
        Mvec e012346() const {return this->extractOneComponent(5,252, 66);}
        Mvec e01234i2() const {return this->extractOneComponent(5,252, 67);}
        Mvec e012356() const {return this->extractOneComponent(5,252, 68);}
        Mvec e01235i2() const {return this->extractOneComponent(5,252, 69);}
        Mvec e01236i2() const {return this->extractOneComponent(5,252, 70);}
        Mvec e012i1024() const {return this->extractOneComponent(5,252, 71);}
        Mvec e012i1025() const {return this->extractOneComponent(5,252, 72);}
        Mvec e012i1026() const {return this->extractOneComponent(5,252, 73);}
        Mvec e012i102i2() const {return this->extractOneComponent(5,252, 74);}
        Mvec e012i145() const {return this->extractOneComponent(5,252, 75);}
        Mvec e012i146() const {return this->extractOneComponent(5,252, 76);}
        Mvec e012i14i2() const {return this->extractOneComponent(5,252, 77);}
        Mvec e012i156() const {return this->extractOneComponent(5,252, 78);}
        Mvec e012i15i2() const {return this->extractOneComponent(5,252, 79);}
        Mvec e012i16i2() const {return this->extractOneComponent(5,252, 80);}
        Mvec e0120245() const {return this->extractOneComponent(5,252, 81);}
        Mvec e0120246() const {return this->extractOneComponent(5,252, 82);}
        Mvec e012024i2() const {return this->extractOneComponent(5,252, 83);}
        Mvec e0120256() const {return this->extractOneComponent(5,252, 84);}
        Mvec e012025i2() const {return this->extractOneComponent(5,252, 85);}
        Mvec e012026i2() const {return this->extractOneComponent(5,252, 86);}
        Mvec e012456() const {return this->extractOneComponent(5,252, 87);}
        Mvec e01245i2() const {return this->extractOneComponent(5,252, 88);}
        Mvec e01246i2() const {return this->extractOneComponent(5,252, 89);}
        Mvec e01256i2() const {return this->extractOneComponent(5,252, 90);}
        Mvec e013i1024() const {return this->extractOneComponent(5,252, 91);}
        Mvec e013i1025() const {return this->extractOneComponent(5,252, 92);}
        Mvec e013i1026() const {return this->extractOneComponent(5,252, 93);}
        Mvec e013i102i2() const {return this->extractOneComponent(5,252, 94);}
        Mvec e013i145() const {return this->extractOneComponent(5,252, 95);}
        Mvec e013i146() const {return this->extractOneComponent(5,252, 96);}
        Mvec e013i14i2() const {return this->extractOneComponent(5,252, 97);}
        Mvec e013i156() const {return this->extractOneComponent(5,252, 98);}
        Mvec e013i15i2() const {return this->extractOneComponent(5,252, 99);}
        Mvec e013i16i2() const {return this->extractOneComponent(5,252, 100);}
        Mvec e0130245() const {return this->extractOneComponent(5,252, 101);}
        Mvec e0130246() const {return this->extractOneComponent(5,252, 102);}
        Mvec e013024i2() const {return this->extractOneComponent(5,252, 103);}
        Mvec e0130256() const {return this->extractOneComponent(5,252, 104);}
        Mvec e013025i2() const {return this->extractOneComponent(5,252, 105);}
        Mvec e013026i2() const {return this->extractOneComponent(5,252, 106);}
        Mvec e013456() const {return this->extractOneComponent(5,252, 107);}
        Mvec e01345i2() const {return this->extractOneComponent(5,252, 108);}
        Mvec e01346i2() const {return this->extractOneComponent(5,252, 109);}
        Mvec e01356i2() const {return this->extractOneComponent(5,252, 110);}
        Mvec e01i10245() const {return this->extractOneComponent(5,252, 111);}
        Mvec e01i10246() const {return this->extractOneComponent(5,252, 112);}
        Mvec e01i1024i2() const {return this->extractOneComponent(5,252, 113);}
        Mvec e01i10256() const {return this->extractOneComponent(5,252, 114);}
        Mvec e01i1025i2() const {return this->extractOneComponent(5,252, 115);}
        Mvec e01i1026i2() const {return this->extractOneComponent(5,252, 116);}
        Mvec e01i1456() const {return this->extractOneComponent(5,252, 117);}
        Mvec e01i145i2() const {return this->extractOneComponent(5,252, 118);}
        Mvec e01i146i2() const {return this->extractOneComponent(5,252, 119);}
        Mvec e01i156i2() const {return this->extractOneComponent(5,252, 120);}
        Mvec e0102456() const {return this->extractOneComponent(5,252, 121);}
        Mvec e010245i2() const {return this->extractOneComponent(5,252, 122);}
        Mvec e010246i2() const {return this->extractOneComponent(5,252, 123);}
        Mvec e010256i2() const {return this->extractOneComponent(5,252, 124);}
        Mvec e01456i2() const {return this->extractOneComponent(5,252, 125);}
        Mvec e123i102() const {return this->extractOneComponent(5,252, 126);}
        Mvec e123i14() const {return this->extractOneComponent(5,252, 127);}
        Mvec e123i15() const {return this->extractOneComponent(5,252, 128);}
        Mvec e123i16() const {return this->extractOneComponent(5,252, 129);}
        Mvec e123i1i2() const {return this->extractOneComponent(5,252, 130);}
        Mvec e123024() const {return this->extractOneComponent(5,252, 131);}
        Mvec e123025() const {return this->extractOneComponent(5,252, 132);}
        Mvec e123026() const {return this->extractOneComponent(5,252, 133);}
        Mvec e12302i2() const {return this->extractOneComponent(5,252, 134);}
        Mvec e12345() const {return this->extractOneComponent(5,252, 135);}
        Mvec e12346() const {return this->extractOneComponent(5,252, 136);}
        Mvec e1234i2() const {return this->extractOneComponent(5,252, 137);}
        Mvec e12356() const {return this->extractOneComponent(5,252, 138);}
        Mvec e1235i2() const {return this->extractOneComponent(5,252, 139);}
        Mvec e1236i2() const {return this->extractOneComponent(5,252, 140);}
        Mvec e12i1024() const {return this->extractOneComponent(5,252, 141);}
        Mvec e12i1025() const {return this->extractOneComponent(5,252, 142);}
        Mvec e12i1026() const {return this->extractOneComponent(5,252, 143);}
        Mvec e12i102i2() const {return this->extractOneComponent(5,252, 144);}
        Mvec e12i145() const {return this->extractOneComponent(5,252, 145);}
        Mvec e12i146() const {return this->extractOneComponent(5,252, 146);}
        Mvec e12i14i2() const {return this->extractOneComponent(5,252, 147);}
        Mvec e12i156() const {return this->extractOneComponent(5,252, 148);}
        Mvec e12i15i2() const {return this->extractOneComponent(5,252, 149);}
        Mvec e12i16i2() const {return this->extractOneComponent(5,252, 150);}
        Mvec e120245() const {return this->extractOneComponent(5,252, 151);}
        Mvec e120246() const {return this->extractOneComponent(5,252, 152);}
        Mvec e12024i2() const {return this->extractOneComponent(5,252, 153);}
        Mvec e120256() const {return this->extractOneComponent(5,252, 154);}
        Mvec e12025i2() const {return this->extractOneComponent(5,252, 155);}
        Mvec e12026i2() const {return this->extractOneComponent(5,252, 156);}
        Mvec e12456() const {return this->extractOneComponent(5,252, 157);}
        Mvec e1245i2() const {return this->extractOneComponent(5,252, 158);}
        Mvec e1246i2() const {return this->extractOneComponent(5,252, 159);}
        Mvec e1256i2() const {return this->extractOneComponent(5,252, 160);}
        Mvec e13i1024() const {return this->extractOneComponent(5,252, 161);}
        Mvec e13i1025() const {return this->extractOneComponent(5,252, 162);}
        Mvec e13i1026() const {return this->extractOneComponent(5,252, 163);}
        Mvec e13i102i2() const {return this->extractOneComponent(5,252, 164);}
        Mvec e13i145() const {return this->extractOneComponent(5,252, 165);}
        Mvec e13i146() const {return this->extractOneComponent(5,252, 166);}
        Mvec e13i14i2() const {return this->extractOneComponent(5,252, 167);}
        Mvec e13i156() const {return this->extractOneComponent(5,252, 168);}
        Mvec e13i15i2() const {return this->extractOneComponent(5,252, 169);}
        Mvec e13i16i2() const {return this->extractOneComponent(5,252, 170);}
        Mvec e130245() const {return this->extractOneComponent(5,252, 171);}
        Mvec e130246() const {return this->extractOneComponent(5,252, 172);}
        Mvec e13024i2() const {return this->extractOneComponent(5,252, 173);}
        Mvec e130256() const {return this->extractOneComponent(5,252, 174);}
        Mvec e13025i2() const {return this->extractOneComponent(5,252, 175);}
        Mvec e13026i2() const {return this->extractOneComponent(5,252, 176);}
        Mvec e13456() const {return this->extractOneComponent(5,252, 177);}
        Mvec e1345i2() const {return this->extractOneComponent(5,252, 178);}
        Mvec e1346i2() const {return this->extractOneComponent(5,252, 179);}
        Mvec e1356i2() const {return this->extractOneComponent(5,252, 180);}
        Mvec e1i10245() const {return this->extractOneComponent(5,252, 181);}
        Mvec e1i10246() const {return this->extractOneComponent(5,252, 182);}
        Mvec e1i1024i2() const {return this->extractOneComponent(5,252, 183);}
        Mvec e1i10256() const {return this->extractOneComponent(5,252, 184);}
        Mvec e1i1025i2() const {return this->extractOneComponent(5,252, 185);}
        Mvec e1i1026i2() const {return this->extractOneComponent(5,252, 186);}
        Mvec e1i1456() const {return this->extractOneComponent(5,252, 187);}
        Mvec e1i145i2() const {return this->extractOneComponent(5,252, 188);}
        Mvec e1i146i2() const {return this->extractOneComponent(5,252, 189);}
        Mvec e1i156i2() const {return this->extractOneComponent(5,252, 190);}
        Mvec e102456() const {return this->extractOneComponent(5,252, 191);}
        Mvec e10245i2() const {return this->extractOneComponent(5,252, 192);}
        Mvec e10246i2() const {return this->extractOneComponent(5,252, 193);}
        Mvec e10256i2() const {return this->extractOneComponent(5,252, 194);}
        Mvec e1456i2() const {return this->extractOneComponent(5,252, 195);}
        Mvec e23i1024() const {return this->extractOneComponent(5,252, 196);}
        Mvec e23i1025() const {return this->extractOneComponent(5,252, 197);}
        Mvec e23i1026() const {return this->extractOneComponent(5,252, 198);}
        Mvec e23i102i2() const {return this->extractOneComponent(5,252, 199);}
        Mvec e23i145() const {return this->extractOneComponent(5,252, 200);}
        Mvec e23i146() const {return this->extractOneComponent(5,252, 201);}
        Mvec e23i14i2() const {return this->extractOneComponent(5,252, 202);}
        Mvec e23i156() const {return this->extractOneComponent(5,252, 203);}
        Mvec e23i15i2() const {return this->extractOneComponent(5,252, 204);}
        Mvec e23i16i2() const {return this->extractOneComponent(5,252, 205);}
        Mvec e230245() const {return this->extractOneComponent(5,252, 206);}
        Mvec e230246() const {return this->extractOneComponent(5,252, 207);}
        Mvec e23024i2() const {return this->extractOneComponent(5,252, 208);}
        Mvec e230256() const {return this->extractOneComponent(5,252, 209);}
        Mvec e23025i2() const {return this->extractOneComponent(5,252, 210);}
        Mvec e23026i2() const {return this->extractOneComponent(5,252, 211);}
        Mvec e23456() const {return this->extractOneComponent(5,252, 212);}
        Mvec e2345i2() const {return this->extractOneComponent(5,252, 213);}
        Mvec e2346i2() const {return this->extractOneComponent(5,252, 214);}
        Mvec e2356i2() const {return this->extractOneComponent(5,252, 215);}
        Mvec e2i10245() const {return this->extractOneComponent(5,252, 216);}
        Mvec e2i10246() const {return this->extractOneComponent(5,252, 217);}
        Mvec e2i1024i2() const {return this->extractOneComponent(5,252, 218);}
        Mvec e2i10256() const {return this->extractOneComponent(5,252, 219);}
        Mvec e2i1025i2() const {return this->extractOneComponent(5,252, 220);}
        Mvec e2i1026i2() const {return this->extractOneComponent(5,252, 221);}
        Mvec e2i1456() const {return this->extractOneComponent(5,252, 222);}
        Mvec e2i145i2() const {return this->extractOneComponent(5,252, 223);}
        Mvec e2i146i2() const {return this->extractOneComponent(5,252, 224);}
        Mvec e2i156i2() const {return this->extractOneComponent(5,252, 225);}
        Mvec e202456() const {return this->extractOneComponent(5,252, 226);}
        Mvec e20245i2() const {return this->extractOneComponent(5,252, 227);}
        Mvec e20246i2() const {return this->extractOneComponent(5,252, 228);}
        Mvec e20256i2() const {return this->extractOneComponent(5,252, 229);}
        Mvec e2456i2() const {return this->extractOneComponent(5,252, 230);}
        Mvec e3i10245() const {return this->extractOneComponent(5,252, 231);}
        Mvec e3i10246() const {return this->extractOneComponent(5,252, 232);}
        Mvec e3i1024i2() const {return this->extractOneComponent(5,252, 233);}
        Mvec e3i10256() const {return this->extractOneComponent(5,252, 234);}
        Mvec e3i1025i2() const {return this->extractOneComponent(5,252, 235);}
        Mvec e3i1026i2() const {return this->extractOneComponent(5,252, 236);}
        Mvec e3i1456() const {return this->extractOneComponent(5,252, 237);}
        Mvec e3i145i2() const {return this->extractOneComponent(5,252, 238);}
        Mvec e3i146i2() const {return this->extractOneComponent(5,252, 239);}
        Mvec e3i156i2() const {return this->extractOneComponent(5,252, 240);}
        Mvec e302456() const {return this->extractOneComponent(5,252, 241);}
        Mvec e30245i2() const {return this->extractOneComponent(5,252, 242);}
        Mvec e30246i2() const {return this->extractOneComponent(5,252, 243);}
        Mvec e30256i2() const {return this->extractOneComponent(5,252, 244);}
        Mvec e3456i2() const {return this->extractOneComponent(5,252, 245);}
        Mvec ei102456() const {return this->extractOneComponent(5,252, 246);}
        Mvec ei10245i2() const {return this->extractOneComponent(5,252, 247);}
        Mvec ei10246i2() const {return this->extractOneComponent(5,252, 248);}
        Mvec ei10256i2() const {return this->extractOneComponent(5,252, 249);}
        Mvec ei1456i2() const {return this->extractOneComponent(5,252, 250);}
        Mvec e02456i2() const {return this->extractOneComponent(5,252, 251);}
        Mvec e01123i102() const {return this->extractOneComponent(6,210, 0);}
        Mvec e01123i14() const {return this->extractOneComponent(6,210, 1);}
        Mvec e01123i15() const {return this->extractOneComponent(6,210, 2);}
        Mvec e01123i16() const {return this->extractOneComponent(6,210, 3);}
        Mvec e01123i1i2() const {return this->extractOneComponent(6,210, 4);}
        Mvec e01123024() const {return this->extractOneComponent(6,210, 5);}
        Mvec e01123025() const {return this->extractOneComponent(6,210, 6);}
        Mvec e01123026() const {return this->extractOneComponent(6,210, 7);}
        Mvec e0112302i2() const {return this->extractOneComponent(6,210, 8);}
        Mvec e0112345() const {return this->extractOneComponent(6,210, 9);}
        Mvec e0112346() const {return this->extractOneComponent(6,210, 10);}
        Mvec e011234i2() const {return this->extractOneComponent(6,210, 11);}
        Mvec e0112356() const {return this->extractOneComponent(6,210, 12);}
        Mvec e011235i2() const {return this->extractOneComponent(6,210, 13);}
        Mvec e011236i2() const {return this->extractOneComponent(6,210, 14);}
        Mvec e0112i1024() const {return this->extractOneComponent(6,210, 15);}
        Mvec e0112i1025() const {return this->extractOneComponent(6,210, 16);}
        Mvec e0112i1026() const {return this->extractOneComponent(6,210, 17);}
        Mvec e0112i102i2() const {return this->extractOneComponent(6,210, 18);}
        Mvec e0112i145() const {return this->extractOneComponent(6,210, 19);}
        Mvec e0112i146() const {return this->extractOneComponent(6,210, 20);}
        Mvec e0112i14i2() const {return this->extractOneComponent(6,210, 21);}
        Mvec e0112i156() const {return this->extractOneComponent(6,210, 22);}
        Mvec e0112i15i2() const {return this->extractOneComponent(6,210, 23);}
        Mvec e0112i16i2() const {return this->extractOneComponent(6,210, 24);}
        Mvec e01120245() const {return this->extractOneComponent(6,210, 25);}
        Mvec e01120246() const {return this->extractOneComponent(6,210, 26);}
        Mvec e0112024i2() const {return this->extractOneComponent(6,210, 27);}
        Mvec e01120256() const {return this->extractOneComponent(6,210, 28);}
        Mvec e0112025i2() const {return this->extractOneComponent(6,210, 29);}
        Mvec e0112026i2() const {return this->extractOneComponent(6,210, 30);}
        Mvec e0112456() const {return this->extractOneComponent(6,210, 31);}
        Mvec e011245i2() const {return this->extractOneComponent(6,210, 32);}
        Mvec e011246i2() const {return this->extractOneComponent(6,210, 33);}
        Mvec e011256i2() const {return this->extractOneComponent(6,210, 34);}
        Mvec e0113i1024() const {return this->extractOneComponent(6,210, 35);}
        Mvec e0113i1025() const {return this->extractOneComponent(6,210, 36);}
        Mvec e0113i1026() const {return this->extractOneComponent(6,210, 37);}
        Mvec e0113i102i2() const {return this->extractOneComponent(6,210, 38);}
        Mvec e0113i145() const {return this->extractOneComponent(6,210, 39);}
        Mvec e0113i146() const {return this->extractOneComponent(6,210, 40);}
        Mvec e0113i14i2() const {return this->extractOneComponent(6,210, 41);}
        Mvec e0113i156() const {return this->extractOneComponent(6,210, 42);}
        Mvec e0113i15i2() const {return this->extractOneComponent(6,210, 43);}
        Mvec e0113i16i2() const {return this->extractOneComponent(6,210, 44);}
        Mvec e01130245() const {return this->extractOneComponent(6,210, 45);}
        Mvec e01130246() const {return this->extractOneComponent(6,210, 46);}
        Mvec e0113024i2() const {return this->extractOneComponent(6,210, 47);}
        Mvec e01130256() const {return this->extractOneComponent(6,210, 48);}
        Mvec e0113025i2() const {return this->extractOneComponent(6,210, 49);}
        Mvec e0113026i2() const {return this->extractOneComponent(6,210, 50);}
        Mvec e0113456() const {return this->extractOneComponent(6,210, 51);}
        Mvec e011345i2() const {return this->extractOneComponent(6,210, 52);}
        Mvec e011346i2() const {return this->extractOneComponent(6,210, 53);}
        Mvec e011356i2() const {return this->extractOneComponent(6,210, 54);}
        Mvec e011i10245() const {return this->extractOneComponent(6,210, 55);}
        Mvec e011i10246() const {return this->extractOneComponent(6,210, 56);}
        Mvec e011i1024i2() const {return this->extractOneComponent(6,210, 57);}
        Mvec e011i10256() const {return this->extractOneComponent(6,210, 58);}
        Mvec e011i1025i2() const {return this->extractOneComponent(6,210, 59);}
        Mvec e011i1026i2() const {return this->extractOneComponent(6,210, 60);}
        Mvec e011i1456() const {return this->extractOneComponent(6,210, 61);}
        Mvec e011i145i2() const {return this->extractOneComponent(6,210, 62);}
        Mvec e011i146i2() const {return this->extractOneComponent(6,210, 63);}
        Mvec e011i156i2() const {return this->extractOneComponent(6,210, 64);}
        Mvec e01102456() const {return this->extractOneComponent(6,210, 65);}
        Mvec e0110245i2() const {return this->extractOneComponent(6,210, 66);}
        Mvec e0110246i2() const {return this->extractOneComponent(6,210, 67);}
        Mvec e0110256i2() const {return this->extractOneComponent(6,210, 68);}
        Mvec e011456i2() const {return this->extractOneComponent(6,210, 69);}
        Mvec e0123i1024() const {return this->extractOneComponent(6,210, 70);}
        Mvec e0123i1025() const {return this->extractOneComponent(6,210, 71);}
        Mvec e0123i1026() const {return this->extractOneComponent(6,210, 72);}
        Mvec e0123i102i2() const {return this->extractOneComponent(6,210, 73);}
        Mvec e0123i145() const {return this->extractOneComponent(6,210, 74);}
        Mvec e0123i146() const {return this->extractOneComponent(6,210, 75);}
        Mvec e0123i14i2() const {return this->extractOneComponent(6,210, 76);}
        Mvec e0123i156() const {return this->extractOneComponent(6,210, 77);}
        Mvec e0123i15i2() const {return this->extractOneComponent(6,210, 78);}
        Mvec e0123i16i2() const {return this->extractOneComponent(6,210, 79);}
        Mvec e01230245() const {return this->extractOneComponent(6,210, 80);}
        Mvec e01230246() const {return this->extractOneComponent(6,210, 81);}
        Mvec e0123024i2() const {return this->extractOneComponent(6,210, 82);}
        Mvec e01230256() const {return this->extractOneComponent(6,210, 83);}
        Mvec e0123025i2() const {return this->extractOneComponent(6,210, 84);}
        Mvec e0123026i2() const {return this->extractOneComponent(6,210, 85);}
        Mvec e0123456() const {return this->extractOneComponent(6,210, 86);}
        Mvec e012345i2() const {return this->extractOneComponent(6,210, 87);}
        Mvec e012346i2() const {return this->extractOneComponent(6,210, 88);}
        Mvec e012356i2() const {return this->extractOneComponent(6,210, 89);}
        Mvec e012i10245() const {return this->extractOneComponent(6,210, 90);}
        Mvec e012i10246() const {return this->extractOneComponent(6,210, 91);}
        Mvec e012i1024i2() const {return this->extractOneComponent(6,210, 92);}
        Mvec e012i10256() const {return this->extractOneComponent(6,210, 93);}
        Mvec e012i1025i2() const {return this->extractOneComponent(6,210, 94);}
        Mvec e012i1026i2() const {return this->extractOneComponent(6,210, 95);}
        Mvec e012i1456() const {return this->extractOneComponent(6,210, 96);}
        Mvec e012i145i2() const {return this->extractOneComponent(6,210, 97);}
        Mvec e012i146i2() const {return this->extractOneComponent(6,210, 98);}
        Mvec e012i156i2() const {return this->extractOneComponent(6,210, 99);}
        Mvec e01202456() const {return this->extractOneComponent(6,210, 100);}
        Mvec e0120245i2() const {return this->extractOneComponent(6,210, 101);}
        Mvec e0120246i2() const {return this->extractOneComponent(6,210, 102);}
        Mvec e0120256i2() const {return this->extractOneComponent(6,210, 103);}
        Mvec e012456i2() const {return this->extractOneComponent(6,210, 104);}
        Mvec e013i10245() const {return this->extractOneComponent(6,210, 105);}
        Mvec e013i10246() const {return this->extractOneComponent(6,210, 106);}
        Mvec e013i1024i2() const {return this->extractOneComponent(6,210, 107);}
        Mvec e013i10256() const {return this->extractOneComponent(6,210, 108);}
        Mvec e013i1025i2() const {return this->extractOneComponent(6,210, 109);}
        Mvec e013i1026i2() const {return this->extractOneComponent(6,210, 110);}
        Mvec e013i1456() const {return this->extractOneComponent(6,210, 111);}
        Mvec e013i145i2() const {return this->extractOneComponent(6,210, 112);}
        Mvec e013i146i2() const {return this->extractOneComponent(6,210, 113);}
        Mvec e013i156i2() const {return this->extractOneComponent(6,210, 114);}
        Mvec e01302456() const {return this->extractOneComponent(6,210, 115);}
        Mvec e0130245i2() const {return this->extractOneComponent(6,210, 116);}
        Mvec e0130246i2() const {return this->extractOneComponent(6,210, 117);}
        Mvec e0130256i2() const {return this->extractOneComponent(6,210, 118);}
        Mvec e013456i2() const {return this->extractOneComponent(6,210, 119);}
        Mvec e01i102456() const {return this->extractOneComponent(6,210, 120);}
        Mvec e01i10245i2() const {return this->extractOneComponent(6,210, 121);}
        Mvec e01i10246i2() const {return this->extractOneComponent(6,210, 122);}
        Mvec e01i10256i2() const {return this->extractOneComponent(6,210, 123);}
        Mvec e01i1456i2() const {return this->extractOneComponent(6,210, 124);}
        Mvec e0102456i2() const {return this->extractOneComponent(6,210, 125);}
        Mvec e123i1024() const {return this->extractOneComponent(6,210, 126);}
        Mvec e123i1025() const {return this->extractOneComponent(6,210, 127);}
        Mvec e123i1026() const {return this->extractOneComponent(6,210, 128);}
        Mvec e123i102i2() const {return this->extractOneComponent(6,210, 129);}
        Mvec e123i145() const {return this->extractOneComponent(6,210, 130);}
        Mvec e123i146() const {return this->extractOneComponent(6,210, 131);}
        Mvec e123i14i2() const {return this->extractOneComponent(6,210, 132);}
        Mvec e123i156() const {return this->extractOneComponent(6,210, 133);}
        Mvec e123i15i2() const {return this->extractOneComponent(6,210, 134);}
        Mvec e123i16i2() const {return this->extractOneComponent(6,210, 135);}
        Mvec e1230245() const {return this->extractOneComponent(6,210, 136);}
        Mvec e1230246() const {return this->extractOneComponent(6,210, 137);}
        Mvec e123024i2() const {return this->extractOneComponent(6,210, 138);}
        Mvec e1230256() const {return this->extractOneComponent(6,210, 139);}
        Mvec e123025i2() const {return this->extractOneComponent(6,210, 140);}
        Mvec e123026i2() const {return this->extractOneComponent(6,210, 141);}
        Mvec e123456() const {return this->extractOneComponent(6,210, 142);}
        Mvec e12345i2() const {return this->extractOneComponent(6,210, 143);}
        Mvec e12346i2() const {return this->extractOneComponent(6,210, 144);}
        Mvec e12356i2() const {return this->extractOneComponent(6,210, 145);}
        Mvec e12i10245() const {return this->extractOneComponent(6,210, 146);}
        Mvec e12i10246() const {return this->extractOneComponent(6,210, 147);}
        Mvec e12i1024i2() const {return this->extractOneComponent(6,210, 148);}
        Mvec e12i10256() const {return this->extractOneComponent(6,210, 149);}
        Mvec e12i1025i2() const {return this->extractOneComponent(6,210, 150);}
        Mvec e12i1026i2() const {return this->extractOneComponent(6,210, 151);}
        Mvec e12i1456() const {return this->extractOneComponent(6,210, 152);}
        Mvec e12i145i2() const {return this->extractOneComponent(6,210, 153);}
        Mvec e12i146i2() const {return this->extractOneComponent(6,210, 154);}
        Mvec e12i156i2() const {return this->extractOneComponent(6,210, 155);}
        Mvec e1202456() const {return this->extractOneComponent(6,210, 156);}
        Mvec e120245i2() const {return this->extractOneComponent(6,210, 157);}
        Mvec e120246i2() const {return this->extractOneComponent(6,210, 158);}
        Mvec e120256i2() const {return this->extractOneComponent(6,210, 159);}
        Mvec e12456i2() const {return this->extractOneComponent(6,210, 160);}
        Mvec e13i10245() const {return this->extractOneComponent(6,210, 161);}
        Mvec e13i10246() const {return this->extractOneComponent(6,210, 162);}
        Mvec e13i1024i2() const {return this->extractOneComponent(6,210, 163);}
        Mvec e13i10256() const {return this->extractOneComponent(6,210, 164);}
        Mvec e13i1025i2() const {return this->extractOneComponent(6,210, 165);}
        Mvec e13i1026i2() const {return this->extractOneComponent(6,210, 166);}
        Mvec e13i1456() const {return this->extractOneComponent(6,210, 167);}
        Mvec e13i145i2() const {return this->extractOneComponent(6,210, 168);}
        Mvec e13i146i2() const {return this->extractOneComponent(6,210, 169);}
        Mvec e13i156i2() const {return this->extractOneComponent(6,210, 170);}
        Mvec e1302456() const {return this->extractOneComponent(6,210, 171);}
        Mvec e130245i2() const {return this->extractOneComponent(6,210, 172);}
        Mvec e130246i2() const {return this->extractOneComponent(6,210, 173);}
        Mvec e130256i2() const {return this->extractOneComponent(6,210, 174);}
        Mvec e13456i2() const {return this->extractOneComponent(6,210, 175);}
        Mvec e1i102456() const {return this->extractOneComponent(6,210, 176);}
        Mvec e1i10245i2() const {return this->extractOneComponent(6,210, 177);}
        Mvec e1i10246i2() const {return this->extractOneComponent(6,210, 178);}
        Mvec e1i10256i2() const {return this->extractOneComponent(6,210, 179);}
        Mvec e1i1456i2() const {return this->extractOneComponent(6,210, 180);}
        Mvec e102456i2() const {return this->extractOneComponent(6,210, 181);}
        Mvec e23i10245() const {return this->extractOneComponent(6,210, 182);}
        Mvec e23i10246() const {return this->extractOneComponent(6,210, 183);}
        Mvec e23i1024i2() const {return this->extractOneComponent(6,210, 184);}
        Mvec e23i10256() const {return this->extractOneComponent(6,210, 185);}
        Mvec e23i1025i2() const {return this->extractOneComponent(6,210, 186);}
        Mvec e23i1026i2() const {return this->extractOneComponent(6,210, 187);}
        Mvec e23i1456() const {return this->extractOneComponent(6,210, 188);}
        Mvec e23i145i2() const {return this->extractOneComponent(6,210, 189);}
        Mvec e23i146i2() const {return this->extractOneComponent(6,210, 190);}
        Mvec e23i156i2() const {return this->extractOneComponent(6,210, 191);}
        Mvec e2302456() const {return this->extractOneComponent(6,210, 192);}
        Mvec e230245i2() const {return this->extractOneComponent(6,210, 193);}
        Mvec e230246i2() const {return this->extractOneComponent(6,210, 194);}
        Mvec e230256i2() const {return this->extractOneComponent(6,210, 195);}
        Mvec e23456i2() const {return this->extractOneComponent(6,210, 196);}
        Mvec e2i102456() const {return this->extractOneComponent(6,210, 197);}
        Mvec e2i10245i2() const {return this->extractOneComponent(6,210, 198);}
        Mvec e2i10246i2() const {return this->extractOneComponent(6,210, 199);}
        Mvec e2i10256i2() const {return this->extractOneComponent(6,210, 200);}
        Mvec e2i1456i2() const {return this->extractOneComponent(6,210, 201);}
        Mvec e202456i2() const {return this->extractOneComponent(6,210, 202);}
        Mvec e3i102456() const {return this->extractOneComponent(6,210, 203);}
        Mvec e3i10245i2() const {return this->extractOneComponent(6,210, 204);}
        Mvec e3i10246i2() const {return this->extractOneComponent(6,210, 205);}
        Mvec e3i10256i2() const {return this->extractOneComponent(6,210, 206);}
        Mvec e3i1456i2() const {return this->extractOneComponent(6,210, 207);}
        Mvec e302456i2() const {return this->extractOneComponent(6,210, 208);}
        Mvec ei102456i2() const {return this->extractOneComponent(6,210, 209);}
        Mvec e01123i1024() const {return this->extractOneComponent(7,120, 0);}
        Mvec e01123i1025() const {return this->extractOneComponent(7,120, 1);}
        Mvec e01123i1026() const {return this->extractOneComponent(7,120, 2);}
        Mvec e01123i102i2() const {return this->extractOneComponent(7,120, 3);}
        Mvec e01123i145() const {return this->extractOneComponent(7,120, 4);}
        Mvec e01123i146() const {return this->extractOneComponent(7,120, 5);}
        Mvec e01123i14i2() const {return this->extractOneComponent(7,120, 6);}
        Mvec e01123i156() const {return this->extractOneComponent(7,120, 7);}
        Mvec e01123i15i2() const {return this->extractOneComponent(7,120, 8);}
        Mvec e01123i16i2() const {return this->extractOneComponent(7,120, 9);}
        Mvec e011230245() const {return this->extractOneComponent(7,120, 10);}
        Mvec e011230246() const {return this->extractOneComponent(7,120, 11);}
        Mvec e01123024i2() const {return this->extractOneComponent(7,120, 12);}
        Mvec e011230256() const {return this->extractOneComponent(7,120, 13);}
        Mvec e01123025i2() const {return this->extractOneComponent(7,120, 14);}
        Mvec e01123026i2() const {return this->extractOneComponent(7,120, 15);}
        Mvec e01123456() const {return this->extractOneComponent(7,120, 16);}
        Mvec e0112345i2() const {return this->extractOneComponent(7,120, 17);}
        Mvec e0112346i2() const {return this->extractOneComponent(7,120, 18);}
        Mvec e0112356i2() const {return this->extractOneComponent(7,120, 19);}
        Mvec e0112i10245() const {return this->extractOneComponent(7,120, 20);}
        Mvec e0112i10246() const {return this->extractOneComponent(7,120, 21);}
        Mvec e0112i1024i2() const {return this->extractOneComponent(7,120, 22);}
        Mvec e0112i10256() const {return this->extractOneComponent(7,120, 23);}
        Mvec e0112i1025i2() const {return this->extractOneComponent(7,120, 24);}
        Mvec e0112i1026i2() const {return this->extractOneComponent(7,120, 25);}
        Mvec e0112i1456() const {return this->extractOneComponent(7,120, 26);}
        Mvec e0112i145i2() const {return this->extractOneComponent(7,120, 27);}
        Mvec e0112i146i2() const {return this->extractOneComponent(7,120, 28);}
        Mvec e0112i156i2() const {return this->extractOneComponent(7,120, 29);}
        Mvec e011202456() const {return this->extractOneComponent(7,120, 30);}
        Mvec e01120245i2() const {return this->extractOneComponent(7,120, 31);}
        Mvec e01120246i2() const {return this->extractOneComponent(7,120, 32);}
        Mvec e01120256i2() const {return this->extractOneComponent(7,120, 33);}
        Mvec e0112456i2() const {return this->extractOneComponent(7,120, 34);}
        Mvec e0113i10245() const {return this->extractOneComponent(7,120, 35);}
        Mvec e0113i10246() const {return this->extractOneComponent(7,120, 36);}
        Mvec e0113i1024i2() const {return this->extractOneComponent(7,120, 37);}
        Mvec e0113i10256() const {return this->extractOneComponent(7,120, 38);}
        Mvec e0113i1025i2() const {return this->extractOneComponent(7,120, 39);}
        Mvec e0113i1026i2() const {return this->extractOneComponent(7,120, 40);}
        Mvec e0113i1456() const {return this->extractOneComponent(7,120, 41);}
        Mvec e0113i145i2() const {return this->extractOneComponent(7,120, 42);}
        Mvec e0113i146i2() const {return this->extractOneComponent(7,120, 43);}
        Mvec e0113i156i2() const {return this->extractOneComponent(7,120, 44);}
        Mvec e011302456() const {return this->extractOneComponent(7,120, 45);}
        Mvec e01130245i2() const {return this->extractOneComponent(7,120, 46);}
        Mvec e01130246i2() const {return this->extractOneComponent(7,120, 47);}
        Mvec e01130256i2() const {return this->extractOneComponent(7,120, 48);}
        Mvec e0113456i2() const {return this->extractOneComponent(7,120, 49);}
        Mvec e011i102456() const {return this->extractOneComponent(7,120, 50);}
        Mvec e011i10245i2() const {return this->extractOneComponent(7,120, 51);}
        Mvec e011i10246i2() const {return this->extractOneComponent(7,120, 52);}
        Mvec e011i10256i2() const {return this->extractOneComponent(7,120, 53);}
        Mvec e011i1456i2() const {return this->extractOneComponent(7,120, 54);}
        Mvec e01102456i2() const {return this->extractOneComponent(7,120, 55);}
        Mvec e0123i10245() const {return this->extractOneComponent(7,120, 56);}
        Mvec e0123i10246() const {return this->extractOneComponent(7,120, 57);}
        Mvec e0123i1024i2() const {return this->extractOneComponent(7,120, 58);}
        Mvec e0123i10256() const {return this->extractOneComponent(7,120, 59);}
        Mvec e0123i1025i2() const {return this->extractOneComponent(7,120, 60);}
        Mvec e0123i1026i2() const {return this->extractOneComponent(7,120, 61);}
        Mvec e0123i1456() const {return this->extractOneComponent(7,120, 62);}
        Mvec e0123i145i2() const {return this->extractOneComponent(7,120, 63);}
        Mvec e0123i146i2() const {return this->extractOneComponent(7,120, 64);}
        Mvec e0123i156i2() const {return this->extractOneComponent(7,120, 65);}
        Mvec e012302456() const {return this->extractOneComponent(7,120, 66);}
        Mvec e01230245i2() const {return this->extractOneComponent(7,120, 67);}
        Mvec e01230246i2() const {return this->extractOneComponent(7,120, 68);}
        Mvec e01230256i2() const {return this->extractOneComponent(7,120, 69);}
        Mvec e0123456i2() const {return this->extractOneComponent(7,120, 70);}
        Mvec e012i102456() const {return this->extractOneComponent(7,120, 71);}
        Mvec e012i10245i2() const {return this->extractOneComponent(7,120, 72);}
        Mvec e012i10246i2() const {return this->extractOneComponent(7,120, 73);}
        Mvec e012i10256i2() const {return this->extractOneComponent(7,120, 74);}
        Mvec e012i1456i2() const {return this->extractOneComponent(7,120, 75);}
        Mvec e01202456i2() const {return this->extractOneComponent(7,120, 76);}
        Mvec e013i102456() const {return this->extractOneComponent(7,120, 77);}
        Mvec e013i10245i2() const {return this->extractOneComponent(7,120, 78);}
        Mvec e013i10246i2() const {return this->extractOneComponent(7,120, 79);}
        Mvec e013i10256i2() const {return this->extractOneComponent(7,120, 80);}
        Mvec e013i1456i2() const {return this->extractOneComponent(7,120, 81);}
        Mvec e01302456i2() const {return this->extractOneComponent(7,120, 82);}
        Mvec e01i102456i2() const {return this->extractOneComponent(7,120, 83);}
        Mvec e123i10245() const {return this->extractOneComponent(7,120, 84);}
        Mvec e123i10246() const {return this->extractOneComponent(7,120, 85);}
        Mvec e123i1024i2() const {return this->extractOneComponent(7,120, 86);}
        Mvec e123i10256() const {return this->extractOneComponent(7,120, 87);}
        Mvec e123i1025i2() const {return this->extractOneComponent(7,120, 88);}
        Mvec e123i1026i2() const {return this->extractOneComponent(7,120, 89);}
        Mvec e123i1456() const {return this->extractOneComponent(7,120, 90);}
        Mvec e123i145i2() const {return this->extractOneComponent(7,120, 91);}
        Mvec e123i146i2() const {return this->extractOneComponent(7,120, 92);}
        Mvec e123i156i2() const {return this->extractOneComponent(7,120, 93);}
        Mvec e12302456() const {return this->extractOneComponent(7,120, 94);}
        Mvec e1230245i2() const {return this->extractOneComponent(7,120, 95);}
        Mvec e1230246i2() const {return this->extractOneComponent(7,120, 96);}
        Mvec e1230256i2() const {return this->extractOneComponent(7,120, 97);}
        Mvec e123456i2() const {return this->extractOneComponent(7,120, 98);}
        Mvec e12i102456() const {return this->extractOneComponent(7,120, 99);}
        Mvec e12i10245i2() const {return this->extractOneComponent(7,120, 100);}
        Mvec e12i10246i2() const {return this->extractOneComponent(7,120, 101);}
        Mvec e12i10256i2() const {return this->extractOneComponent(7,120, 102);}
        Mvec e12i1456i2() const {return this->extractOneComponent(7,120, 103);}
        Mvec e1202456i2() const {return this->extractOneComponent(7,120, 104);}
        Mvec e13i102456() const {return this->extractOneComponent(7,120, 105);}
        Mvec e13i10245i2() const {return this->extractOneComponent(7,120, 106);}
        Mvec e13i10246i2() const {return this->extractOneComponent(7,120, 107);}
        Mvec e13i10256i2() const {return this->extractOneComponent(7,120, 108);}
        Mvec e13i1456i2() const {return this->extractOneComponent(7,120, 109);}
        Mvec e1302456i2() const {return this->extractOneComponent(7,120, 110);}
        Mvec e1i102456i2() const {return this->extractOneComponent(7,120, 111);}
        Mvec e23i102456() const {return this->extractOneComponent(7,120, 112);}
        Mvec e23i10245i2() const {return this->extractOneComponent(7,120, 113);}
        Mvec e23i10246i2() const {return this->extractOneComponent(7,120, 114);}
        Mvec e23i10256i2() const {return this->extractOneComponent(7,120, 115);}
        Mvec e23i1456i2() const {return this->extractOneComponent(7,120, 116);}
        Mvec e2302456i2() const {return this->extractOneComponent(7,120, 117);}
        Mvec e2i102456i2() const {return this->extractOneComponent(7,120, 118);}
        Mvec e3i102456i2() const {return this->extractOneComponent(7,120, 119);}
        Mvec e01123i10245() const {return this->extractOneComponent(8,45, 0);}
        Mvec e01123i10246() const {return this->extractOneComponent(8,45, 1);}
        Mvec e01123i1024i2() const {return this->extractOneComponent(8,45, 2);}
        Mvec e01123i10256() const {return this->extractOneComponent(8,45, 3);}
        Mvec e01123i1025i2() const {return this->extractOneComponent(8,45, 4);}
        Mvec e01123i1026i2() const {return this->extractOneComponent(8,45, 5);}
        Mvec e01123i1456() const {return this->extractOneComponent(8,45, 6);}
        Mvec e01123i145i2() const {return this->extractOneComponent(8,45, 7);}
        Mvec e01123i146i2() const {return this->extractOneComponent(8,45, 8);}
        Mvec e01123i156i2() const {return this->extractOneComponent(8,45, 9);}
        Mvec e0112302456() const {return this->extractOneComponent(8,45, 10);}
        Mvec e011230245i2() const {return this->extractOneComponent(8,45, 11);}
        Mvec e011230246i2() const {return this->extractOneComponent(8,45, 12);}
        Mvec e011230256i2() const {return this->extractOneComponent(8,45, 13);}
        Mvec e01123456i2() const {return this->extractOneComponent(8,45, 14);}
        Mvec e0112i102456() const {return this->extractOneComponent(8,45, 15);}
        Mvec e0112i10245i2() const {return this->extractOneComponent(8,45, 16);}
        Mvec e0112i10246i2() const {return this->extractOneComponent(8,45, 17);}
        Mvec e0112i10256i2() const {return this->extractOneComponent(8,45, 18);}
        Mvec e0112i1456i2() const {return this->extractOneComponent(8,45, 19);}
        Mvec e011202456i2() const {return this->extractOneComponent(8,45, 20);}
        Mvec e0113i102456() const {return this->extractOneComponent(8,45, 21);}
        Mvec e0113i10245i2() const {return this->extractOneComponent(8,45, 22);}
        Mvec e0113i10246i2() const {return this->extractOneComponent(8,45, 23);}
        Mvec e0113i10256i2() const {return this->extractOneComponent(8,45, 24);}
        Mvec e0113i1456i2() const {return this->extractOneComponent(8,45, 25);}
        Mvec e011302456i2() const {return this->extractOneComponent(8,45, 26);}
        Mvec e011i102456i2() const {return this->extractOneComponent(8,45, 27);}
        Mvec e0123i102456() const {return this->extractOneComponent(8,45, 28);}
        Mvec e0123i10245i2() const {return this->extractOneComponent(8,45, 29);}
        Mvec e0123i10246i2() const {return this->extractOneComponent(8,45, 30);}
        Mvec e0123i10256i2() const {return this->extractOneComponent(8,45, 31);}
        Mvec e0123i1456i2() const {return this->extractOneComponent(8,45, 32);}
        Mvec e012302456i2() const {return this->extractOneComponent(8,45, 33);}
        Mvec e012i102456i2() const {return this->extractOneComponent(8,45, 34);}
        Mvec e013i102456i2() const {return this->extractOneComponent(8,45, 35);}
        Mvec e123i102456() const {return this->extractOneComponent(8,45, 36);}
        Mvec e123i10245i2() const {return this->extractOneComponent(8,45, 37);}
        Mvec e123i10246i2() const {return this->extractOneComponent(8,45, 38);}
        Mvec e123i10256i2() const {return this->extractOneComponent(8,45, 39);}
        Mvec e123i1456i2() const {return this->extractOneComponent(8,45, 40);}
        Mvec e12302456i2() const {return this->extractOneComponent(8,45, 41);}
        Mvec e12i102456i2() const {return this->extractOneComponent(8,45, 42);}
        Mvec e13i102456i2() const {return this->extractOneComponent(8,45, 43);}
        Mvec e23i102456i2() const {return this->extractOneComponent(8,45, 44);}
        Mvec e01123i102456() const {return this->extractOneComponent(9,10, 0);}
        Mvec e01123i10245i2() const {return this->extractOneComponent(9,10, 1);}
        Mvec e01123i10246i2() const {return this->extractOneComponent(9,10, 2);}
        Mvec e01123i10256i2() const {return this->extractOneComponent(9,10, 3);}
        Mvec e01123i1456i2() const {return this->extractOneComponent(9,10, 4);}
        Mvec e0112302456i2() const {return this->extractOneComponent(9,10, 5);}
        Mvec e0112i102456i2() const {return this->extractOneComponent(9,10, 6);}
        Mvec e0113i102456i2() const {return this->extractOneComponent(9,10, 7);}
        Mvec e0123i102456i2() const {return this->extractOneComponent(9,10, 8);}
        Mvec e123i102456i2() const {return this->extractOneComponent(9,10, 9);}
        Mvec e01123i102456i2() const {return this->extractOneComponent(10,1, 0);}


    };  // end of the class definition


    /* ------------------------------------------------------------------------------------------------ */


    template<typename T>
    Mvec<T>::Mvec():gradeBitmap(0)
    {}


    template<typename T>
    Mvec<T>::Mvec(const Mvec& mv) : mvData(mv.mvData), gradeBitmap(mv.gradeBitmap)
    {
        //std::cout << "copy constructor " << std::endl;
    }


    template<typename T>
    Mvec<T>::Mvec(Mvec<T>&& multivector) : mvData(std::move(multivector.mvData)), gradeBitmap(multivector.gradeBitmap)
    {
        // the move constructor
        //std::cout << "move constructor" << std::endl;
    }


    template<typename T>
    template<typename U>
    Mvec<T>::Mvec(const Mvec<U> &mv) : gradeBitmap(mv.gradeBitmap)
    {
        for(auto it=mv.mvData.begin(); it != mv.mvData.end(); ++it){
            Kvec<T> kvec;
            kvec.grade = it->grade;
            kvec.vec = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(binomialArray[it->grade]);
            for(unsigned int i=0; i<it->vec.size(); ++i)
                kvec.vec.coeffRef(i) = T(it->vec.coeff(i));
            mvData.push_back(kvec);
        }
    }


    template<typename T>
    template<typename U>
    Mvec<T>::Mvec(const U val) {
        if(val != U(0)) {
            gradeBitmap = 1;
            Kvec<T> kvec;
            kvec.vec =Eigen::Matrix<T, Eigen::Dynamic, 1>(1);
            kvec.grade =0;
            mvData.push_back(kvec);
            mvData.begin()->vec.coeffRef(0) = val;
        }
    }


    template<typename T>
    Mvec<T>::~Mvec()
    {}


    template<typename T>
    Mvec<T>& Mvec<T>::operator=(const Mvec& mv){
        if(&mv == this) return *this;
        gradeBitmap = mv.gradeBitmap;
        mvData = mv.mvData;
        return *this;
    }


    template<typename T>
    Mvec<T>& Mvec<T>::operator=(Mvec&& mv){
        mvData = std::move(mv.mvData);
        gradeBitmap = mv.gradeBitmap;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator+(const Mvec<T> &mv2) const {
        Mvec<T> mv3(*this);
        for(auto & itMv : mv2.mvData) {
            auto it = mv3.createVectorXdIfDoesNotExist(itMv.grade);
            it->vec += itMv.vec;
        }
        return mv3;
    }


    template<typename U, typename S>
    Mvec<U> operator+(const S &value, const Mvec<U> &mv){
        return mv + value;
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator+(const S &value) const {
        Mvec<T> mv(*this);
        if(value != T(0)) {
            auto it = mv.createVectorXdIfDoesNotExist(0);
            it->vec.coeffRef(0) += value;
        }
        return mv;
    }


    template<typename T>
    Mvec<T>& Mvec<T>::operator+=(const Mvec& mv){
        *this = *this + mv;
        return *this;
    }


    template<typename T>
    Mvec<T> operator-(const Mvec<T> &mv) { // unary -
        Mvec<T> mv2(mv);
        for(auto & itMv : mv2.mvData)
            itMv.vec = -itMv.vec;
        return mv2;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator-(const Mvec<T> &mv2) const {
        Mvec<T> mv3(*this);
        for(auto & itMv : mv2.mvData) {
            auto it = mv3.createVectorXdIfDoesNotExist(itMv.grade);
            it->vec -= itMv.vec;
        }
        return mv3;
    }


    template<typename T>
    Mvec<T>& Mvec<T>::operator-=(const Mvec& mv){
        *this = *this - mv;
        return *this;
    }


    template<typename U, typename S>
    Mvec<U> operator-(const S &value, const Mvec<U> &mv){
        return (-mv) + value;
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator-(const S &value) const {
        Mvec<T> mv(*this);
        if(value != T(0)) {
            auto it = mv.createVectorXdIfDoesNotExist(0);
            it->vec.coeffRef(0) -= value;
        }
        return mv;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator^(const Mvec<T> &mv2) const {
#if 0 // only recursive version
        // Loop over non-empty grade of mv1 and mv2
        // This version (with recursive call) is only faster than the standard recursive call if
        // the multivector mv1 and mv2 are homogeneous or near from homogeneous.
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)    // all per-grade component of mv1
            for(const auto & itMv2 : mv2.mvData)  // all per-grade component of mv2
            {
                unsigned int grade_mv3 = itMv1.grade + itMv2.grade;
                if(grade_mv3 <= algebraDimension) {
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(grade_mv3);
                    outerProductHomogeneous(itMv1.vec, itMv2.vec, itMv3->vec,
                                            itMv1.grade, itMv2.grade, grade_mv3);
                }
            }
        return mv3;
#else // use the adaptative pointer function array
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled outer function using the functions pointer called outerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){
                if(itMv1.grade + itMv2.grade <= (int) algebraDimension ){
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(itMv1.grade + itMv2.grade);
                    outerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                }
            }
        return mv3;
#endif
    }

    template<typename U, typename S>
    Mvec<U> operator^(const S &value, const Mvec<U> &mv) {
        return  mv^value;
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator^(const S &value) const {
        Mvec<T> mv2(*this);
        for(auto & itMv : mv2.mvData)
            itMv.vec *= T(value);
        return mv2;
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator^=(const Mvec &mv) {
        *this = *this ^ mv;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator|(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }


    template<typename U, typename S>
    Mvec<U> operator|(const S &value, const Mvec<U> &mv) {
        return mv | value;
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator|(const S &value) const {
        return (*this) > value ; // inner between a mv and a scalar gives 0 (empty multivector)
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator|=(const Mvec &mv) {
        *this = *this | mv;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator>(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                // right contraction constraint: gradeMv1 >= gradeMv2
                if(itMv1.grade < itMv2.grade)
                    continue;

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }

        return mv3;
    }


    template<typename U, typename S>
    Mvec<U> operator>(const S &value, const Mvec<U> &mv) {
        if( (mv.gradeBitmap & 1) == 0) return Mvec<U>();
        else return Mvec<U>( U(mv.mvData.begin()->vec.coeff(0) * value ) );
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator>(const S &value) const {
        // return mv x value
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            itMv.vec *= T(value);
        return mv;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator<(const Mvec<T> &mv2) const{

        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                // left contraction constraint: gradeMv1 <= gradeMv2
                if(itMv1.grade > itMv2.grade)
                    continue;

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }


    template<typename U, typename S>
    Mvec<U> operator<(const S &value, const Mvec<U> &mv) {
        // return mv x value
        Mvec<U> mv3(mv);
        for(auto & itMv : mv3.mvData)
            itMv.vec *= U(value);
        return mv3;
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator<(const S &value) const {
        if( (gradeBitmap & 1) == 0) return Mvec<T>();
        else return Mvec<T>( T(mvData.begin()->vec.coeff(0) * value ) );
    }


    template<typename T>
    Mvec<T> Mvec<T>::scalarProduct(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                if(itMv1.grade != itMv2.grade)
                    continue;

                // perform the inner product
                int absGradeMv3 = 0;
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }



    template<typename T>
    Mvec<T> Mvec<T>::outerPrimalDual(const Mvec<T> &mv2) const{
        Mvec<T> mv3;

        for(const auto & itMv1 : this->mvData)    // all per-grade component of mv1
            for(const auto & itMv2 : mv2.mvData)  // all per-grade component of mv2
            {
                unsigned int grade_mv3 = itMv1.grade + (algebraDimension-itMv2.grade);
                if(grade_mv3 <= algebraDimension) {
                    // handle the scalar product as well as the left contraction
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(grade_mv3);
                    outerProductPrimalDual(itMv1.vec, itMv2.vec, itMv3->vec,
                                           itMv1.grade, itMv2.grade, (unsigned)(algebraDimension-grade_mv3));

                    // check if the result is non-zero
                    if(!((itMv3->vec.array() != 0.0).any())){
                        mv3.mvData.erase(itMv3);
                        mv3.gradeBitmap &= ~(1<<grade_mv3);
                    }
                }
            }

        return mv3;
    }

    template<typename T>
    Mvec<T> Mvec<T>::outerDualPrimal(const Mvec<T> &mv2) const{
        Mvec<T> mv3;

        for(const auto & itMv1 : this->mvData)    // all per-grade component of mv1
            for(const auto & itMv2 : mv2.mvData)  // all per-grade component of mv2
            {
                unsigned int grade_mv3 = itMv1.grade + (algebraDimension-itMv2.grade);
                if(grade_mv3 <= algebraDimension) {
                    // handle the scalar product as well as the left contraction
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(grade_mv3);
                    outerProductDualPrimal(itMv1.vec, itMv2.vec, itMv3->vec,
                                           itMv1.grade, itMv2.grade, (unsigned)(algebraDimension-grade_mv3));

                    // check if the result is non-zero
                    if(!((itMv3->vec.array() != 0.0).any())){
                        mv3.mvData.erase(itMv3);
                        mv3.gradeBitmap &= ~(1<<grade_mv3);
                    }
                }
            }


        return mv3;

    }

    template<typename T>
    Mvec<T> Mvec<T>::outerDualDual(const Mvec<T> &mv2) const{
        Mvec<T> mv3;

        for(const auto & itMv1 : this->mvData)    // all per-grade component of mv1
            for(const auto & itMv2 : mv2.mvData)  // all per-grade component of mv2
            {
                unsigned int grade_mv3 = itMv1.grade + (algebraDimension-itMv2.grade);
                if(grade_mv3 <= algebraDimension) {
                    // handle the scalar product as well as the left contraction
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(grade_mv3);
                    outerProductDualDual(itMv1.vec, itMv2.vec, itMv3->vec,
                                           itMv1.grade, itMv2.grade, (unsigned)(algebraDimension-grade_mv3));

                    // check if the result is non-zero
                    if(!((itMv3->vec.array() != 0.0).any())){
                        mv3.mvData.erase(itMv3);
                        mv3.gradeBitmap &= ~(1<<grade_mv3);
                    }
                }
            }


        return mv3;

    }



    template<typename T>
    Mvec<T> Mvec<T>::hestenesProduct(const Mvec<T> &mv2) const{
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled inner function using the functions pointer called innerFunctionsContainer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){

                if(itMv1.grade*itMv2.grade == 0)
                    continue;

                // perform the inner product
                int absGradeMv3 = std::abs((int)(itMv1.grade - itMv2.grade));
                auto itMv3 = mv3.createVectorXdIfDoesNotExist(absGradeMv3);
                innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);

                // check if the result is non-zero
                if(!((itMv3->vec.array() != 0.0).any())){
                    mv3.mvData.erase(itMv3);
                    mv3.gradeBitmap &= ~(1<<absGradeMv3);
                }
            }
        return mv3;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator*(const Mvec<T> &mv2) const {
        // Loop over non-empty grade of mv1 and mv2
        // call the right explicit unrolled product function using the functions pointer
        Mvec<T> mv3;
        for(const auto & itMv1 : this->mvData)
            for(const auto & itMv2 : mv2.mvData){
                
                // outer product block
                unsigned int gradeOuter = itMv1.grade + itMv2.grade;
                if(gradeOuter <=  algebraDimension ){
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(gradeOuter);
                    outerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                    if(!((itMv3->vec.array() != 0.0).any())){
                        mv3.mvData.erase(itMv3);
                        mv3.gradeBitmap &= ~(1<<gradeOuter);
                    }
                }

                // inner product block
                unsigned int gradeInner = (unsigned int)std::abs((int)(itMv1.grade-itMv2.grade));
                // when the grade of one of the kvectors is zero, the inner product is the same as the outer product
                if(gradeInner != gradeOuter) {
                    auto itMv3 = mv3.createVectorXdIfDoesNotExist(gradeInner);
                    innerFunctionsContainer<T>[itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                    // check if the result is non-zero
                    if(!((itMv3->vec.array() != 0.0).any())){
                        mv3.mvData.erase(itMv3);
                        mv3.gradeBitmap &= ~(1<<gradeInner);
                    }

                    // geometric product part
                    int gradeMax = std::min(((2*algebraDimension)-gradeOuter)+1,gradeOuter);
                    for (int gradeResult = gradeInner+2; gradeResult < gradeMax; gradeResult+=2) {
                        auto itMv3 = mv3.createVectorXdIfDoesNotExist(gradeResult);
                        geometricFunctionsContainer<T>[gradeResult][itMv1.grade][itMv2.grade](itMv1.vec, itMv2.vec, itMv3->vec);
                        // check if the result is non-zero
                        if(!((itMv3->vec.array() != 0.0).any())){
                            mv3.mvData.erase(itMv3);
                            mv3.gradeBitmap &= ~(1<<gradeResult);
                        }
                    }
                }
            }
        return mv3;
    }


    template<typename U, typename S>
    Mvec<U> operator*(const S &value, const Mvec<U> &mv){
        return mv * value;
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator*(const S &value) const{
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            itMv.vec *= T(value);
        return mv;
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator*=(const Mvec &mv) {
        *this = *this * mv;
        return *this;
    }


    template<typename T>
    Mvec<T> Mvec<T>::operator/(const Mvec<T> &mv2) const {
        return *this * mv2.inv();
    }


    template<typename U, typename S>
    Mvec<U> operator/(const S &value, const Mvec<U> &mv){
        return value * mv.inv();
    }


    template<typename T>
    template<typename S>
    Mvec<T> Mvec<T>::operator/(const S &value) const {
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            itMv.vec /= value;
        return mv;
    }


    template<typename T>
    Mvec<T> &Mvec<T>::operator/=(const Mvec &mv) {
        *this = *this / mv;
        return *this;
    }


    template<typename T>
    Mvec<T> operator~(const Mvec<T> &mv){
        return mv.reverse();
    }


    template<typename T>
    Mvec<T> Mvec<T>::inv() const {
        T n = this->quadraticNorm();
        if(n<std::numeric_limits<T>::epsilon() && n>-std::numeric_limits<T>::epsilon())
            return Mvec<T>(); // return 0, this is was gaviewer does.
        return this->reverse() / n;
    };


    template<typename U>
    Mvec<U> operator!(const Mvec<U> &mv) {
        return mv.dual();
    }

    /// \brief defines the left contraction between two multivectors
    /// \param mv1 - a multivector
    /// \param mv2 - second operand, corresponds to a multivector
    /// \return mv1 left contraction mv2
    template<typename T>
    Mvec<T> leftContraction(const Mvec<T> &mv1, const Mvec<T> &mv2){
        return mv1 < mv2;
    }


    /// \brief defines the right contraction between two multivectors
    /// \param mv1 - a multivector
    /// \param mv2 - second operand, corresponds to a multivector
    /// \return mv1 right contraction mv2
    template<typename T>
    Mvec<T> rightContraction(const Mvec<T> &mv1, const Mvec<T> &mv2){
        return mv1 > mv2;
    }


    /// \brief returns a multivector that only contains the coefficient associated to the pseudoscalar.
    /// \return an empty Mvec if the requested element is not part of the multivector, or the multivector that contains only this element if present in the current multivector.
    template<typename T>
    Mvec<T> I() {
        Mvec<T> mvec;
        mvec[1023] = 1.0;
        return mvec;
    }


    /// \brief return the inverse of the pseudo scalar.
    /// \return returns a multivector corresponding to the inverse of the pseudo scalar.
    template<typename T>
    Mvec<T> Iinv() {
        Mvec<T> mvec;
        mvec[1023] = pseudoScalarInverse; // we do not return just a scalar of type T because the grade should be dimension and not 0.
        return mvec;
    }


    // compute the dual of a multivector (i.e mv* = mv.reverse() * Iinv)
    template<typename T>
    Mvec<T> Mvec<T>::dual() const {
        Mvec<T> mvResult;
        // for each k-vectors of the multivector
        for(auto itMv=mvData.rbegin(); itMv!=mvData.rend(); ++itMv){

            // create the dual k-vector
            Kvec<T> kvec ={itMv->vec,algebraDimension-(itMv->grade)};

            // some elements need to be permuted
            for(unsigned int i=0;i<binomialArray[itMv->grade];++i)
                kvec.vec.coeffRef(dualPermutations[itMv->grade][i]) = itMv->vec.coeff(i);

            // the inner product may involve some constant multiplucation for the dual elements
            kvec.vec = kvec.vec.cwiseProduct(dualCoefficients[itMv->grade].template cast<T>());

            // add the k-vector to the resulting multivector
            mvResult.mvData.push_back(kvec);
            mvResult.gradeBitmap |= (1 << kvec.grade);
        }
        return mvResult;
    };

    // \brief compute the reverse of a multivector
    // \return - the reverse of the multivector
    template<typename T>
    Mvec<T> Mvec<T>::reverse() const {
        Mvec<T> mv(*this);
        for(auto & itMv : mv.mvData)
            if(signReversePerGrade[itMv.grade] == -1)
                itMv.vec *= -1;
        return mv;
    }


    template<typename T>
    Mvec<T> Mvec<T>::grade(const int i) const{

        auto it = findGrade(i);
        Mvec<T> mv;

        // if not found, return the empty multivector
        if(it == mvData.end())
            return mv;

        // else return the grade 'i' data
        mv.mvData.push_back(*it);
        mv.gradeBitmap = 1 << (i);

        return mv;
    }


    template<typename T>
    std::vector<unsigned int> Mvec<T>::grades() const {
        std::vector<unsigned int> mvGrades;
        for(unsigned int i=0; i< algebraDimension+1; ++i){

            if( (gradeBitmap & (1<<i)))
                mvGrades.push_back(i);
        }

        return mvGrades;
    }


    template<typename T>
    void Mvec<T>::roundZero(const T epsilon) {
        // loop over each k-vector of the multivector
        auto itMv=mvData.begin();
        while(itMv != mvData.end()){
            // loop over each element of the k-vector
            for(unsigned int i=0; i<(unsigned int)itMv->vec.size(); ++i)
                if(fabs(itMv->vec.coeff(i)) <= epsilon)
                    itMv->vec.coeffRef(i) = 0.0;
            // if the k-vector is full of 0, remove it
            if(!((itMv->vec.array() != 0.0).any())){
                gradeBitmap = gradeBitmap - (1 << itMv->grade);
                mvData.erase(itMv++);
            }
            else ++itMv;
        }
    }


    template<typename T>
    void Mvec<T>::clear(const int grade) {

        // full erase
        if(grade < 0){
            mvData.clear();
            gradeBitmap = 0;
            return;
        }

        // partial erase
        auto iter = findGrade(grade);
        if(iter != mvData.end()) {
            gradeBitmap = gradeBitmap - (1 << iter->grade);
            iter = mvData.erase(iter);
        }
    }


    template<typename T>
    void Mvec<T>::display() const {
        if(mvData.size() == 0)
            std::cout << " null grade , null value " <<std::endl;

        for(auto itmv=mvData.begin(); itmv!=mvData.end(); ++itmv){
            std::cout << "  grade   : " << itmv->grade  << std::endl;
            std::cout << "  kvector : " << itmv->vec.transpose() << std::endl;
        }
        std::cout << std::endl;
    }


    /// \cond DEV
    template<typename U>
    void traverseKVector(std::ostream &stream, const Eigen::Matrix<U, Eigen::Dynamic, 1> kvector, unsigned int gradeMV, bool& moreThanOne ){

        // version with XOR indices
        // for the current grade, generate all the XOR indices
        std::vector<bool> booleanXorIndex(algebraDimension);
        std::fill(booleanXorIndex.begin(), booleanXorIndex.end(), false);
        // build the first xorIndex of grade 'grade' (with 'grade' values to true).
        std::fill(booleanXorIndex.begin(), booleanXorIndex.begin() + gradeMV, true);
        unsigned int positionInKVector = 0;

        do {
            // convert the vector of bool to the string containing all the basis vectors
            std::string basisBlades={};
            for(unsigned int i=0; i<algebraDimension; ++i) {
                if(booleanXorIndex[i]) {
                    basisBlades += basisVectors[i];
                }
            }

            if(kvector.coeff(positionInKVector)!= 0) {
                if(!(moreThanOne)){
                    stream<< kvector.coeff(positionInKVector) << "*e"+ basisBlades;
                    moreThanOne = true;
                }else{
                    if(kvector.coeff(positionInKVector)>0)
                        stream<< " + " << kvector.coeff(positionInKVector) << "*e" + basisBlades;
                    else
                        stream<< " - " << -kvector.coeff(positionInKVector) << "*e" + basisBlades;
                }
            }

            positionInKVector++;

        } while(std::prev_permutation(booleanXorIndex.begin(), booleanXorIndex.end())); // compute next permutation of the true values of booleanXorIndex

    }
    /// \endcond


    template<typename U>
    std::ostream &operator<<(std::ostream &stream, const Mvec<U> &mvec) {
        if(mvec.mvData.size()==0){
            stream << " 0 ";
            return stream;
        }

        bool moreThanOne = false;
        auto mvIterator = mvec.mvData.begin();

        // if the iterator corresponds to the scalar
        if(mvIterator->grade == 0){
            stream << mvIterator->vec.coeff(0);
            moreThanOne = true;
            ++mvIterator;
        }

        // for all other k-vectors of mvec
        for(; mvIterator!=mvec.mvData.end(); ++mvIterator){

            // call the function that covers a single k-vector
            traverseKVector(stream,mvIterator->vec,mvIterator->grade, moreThanOne);
        }

        if(!moreThanOne)
            stream << " 0 ";

        return stream;
    }

    // static unit multivector functions (returns a multivector with only one component)
    /// \brief return a multivector that contains only the unit basis k-vector 01.
    template<typename T>
    static Mvec<T> e01(){
        Mvec<T> res;
        return res.componentToOne(1, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1.
    template<typename T>
    static Mvec<T> e1(){
        Mvec<T> res;
        return res.componentToOne(1, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2.
    template<typename T>
    static Mvec<T> e2(){
        Mvec<T> res;
        return res.componentToOne(1, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3.
    template<typename T>
    static Mvec<T> e3(){
        Mvec<T> res;
        return res.componentToOne(1, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1.
    template<typename T>
    static Mvec<T> ei1(){
        Mvec<T> res;
        return res.componentToOne(1, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 02.
    template<typename T>
    static Mvec<T> e02(){
        Mvec<T> res;
        return res.componentToOne(1, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 4.
    template<typename T>
    static Mvec<T> e4(){
        Mvec<T> res;
        return res.componentToOne(1, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 5.
    template<typename T>
    static Mvec<T> e5(){
        Mvec<T> res;
        return res.componentToOne(1, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 6.
    template<typename T>
    static Mvec<T> e6(){
        Mvec<T> res;
        return res.componentToOne(1, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i2.
    template<typename T>
    static Mvec<T> ei2(){
        Mvec<T> res;
        return res.componentToOne(1, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011.
    template<typename T>
    static Mvec<T> e011(){
        Mvec<T> res;
        return res.componentToOne(2, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012.
    template<typename T>
    static Mvec<T> e012(){
        Mvec<T> res;
        return res.componentToOne(2, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013.
    template<typename T>
    static Mvec<T> e013(){
        Mvec<T> res;
        return res.componentToOne(2, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1.
    template<typename T>
    static Mvec<T> e01i1(){
        Mvec<T> res;
        return res.componentToOne(2, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0102.
    template<typename T>
    static Mvec<T> e0102(){
        Mvec<T> res;
        return res.componentToOne(2, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 014.
    template<typename T>
    static Mvec<T> e014(){
        Mvec<T> res;
        return res.componentToOne(2, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 015.
    template<typename T>
    static Mvec<T> e015(){
        Mvec<T> res;
        return res.componentToOne(2, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 016.
    template<typename T>
    static Mvec<T> e016(){
        Mvec<T> res;
        return res.componentToOne(2, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i2.
    template<typename T>
    static Mvec<T> e01i2(){
        Mvec<T> res;
        return res.componentToOne(2, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12.
    template<typename T>
    static Mvec<T> e12(){
        Mvec<T> res;
        return res.componentToOne(2, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13.
    template<typename T>
    static Mvec<T> e13(){
        Mvec<T> res;
        return res.componentToOne(2, 10);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1.
    template<typename T>
    static Mvec<T> e1i1(){
        Mvec<T> res;
        return res.componentToOne(2, 11);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 102.
    template<typename T>
    static Mvec<T> e102(){
        Mvec<T> res;
        return res.componentToOne(2, 12);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 14.
    template<typename T>
    static Mvec<T> e14(){
        Mvec<T> res;
        return res.componentToOne(2, 13);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 15.
    template<typename T>
    static Mvec<T> e15(){
        Mvec<T> res;
        return res.componentToOne(2, 14);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 16.
    template<typename T>
    static Mvec<T> e16(){
        Mvec<T> res;
        return res.componentToOne(2, 15);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i2.
    template<typename T>
    static Mvec<T> e1i2(){
        Mvec<T> res;
        return res.componentToOne(2, 16);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23.
    template<typename T>
    static Mvec<T> e23(){
        Mvec<T> res;
        return res.componentToOne(2, 17);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1.
    template<typename T>
    static Mvec<T> e2i1(){
        Mvec<T> res;
        return res.componentToOne(2, 18);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 202.
    template<typename T>
    static Mvec<T> e202(){
        Mvec<T> res;
        return res.componentToOne(2, 19);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 24.
    template<typename T>
    static Mvec<T> e24(){
        Mvec<T> res;
        return res.componentToOne(2, 20);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 25.
    template<typename T>
    static Mvec<T> e25(){
        Mvec<T> res;
        return res.componentToOne(2, 21);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 26.
    template<typename T>
    static Mvec<T> e26(){
        Mvec<T> res;
        return res.componentToOne(2, 22);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i2.
    template<typename T>
    static Mvec<T> e2i2(){
        Mvec<T> res;
        return res.componentToOne(2, 23);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1.
    template<typename T>
    static Mvec<T> e3i1(){
        Mvec<T> res;
        return res.componentToOne(2, 24);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 302.
    template<typename T>
    static Mvec<T> e302(){
        Mvec<T> res;
        return res.componentToOne(2, 25);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 34.
    template<typename T>
    static Mvec<T> e34(){
        Mvec<T> res;
        return res.componentToOne(2, 26);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 35.
    template<typename T>
    static Mvec<T> e35(){
        Mvec<T> res;
        return res.componentToOne(2, 27);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 36.
    template<typename T>
    static Mvec<T> e36(){
        Mvec<T> res;
        return res.componentToOne(2, 28);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i2.
    template<typename T>
    static Mvec<T> e3i2(){
        Mvec<T> res;
        return res.componentToOne(2, 29);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i102.
    template<typename T>
    static Mvec<T> ei102(){
        Mvec<T> res;
        return res.componentToOne(2, 30);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i14.
    template<typename T>
    static Mvec<T> ei14(){
        Mvec<T> res;
        return res.componentToOne(2, 31);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i15.
    template<typename T>
    static Mvec<T> ei15(){
        Mvec<T> res;
        return res.componentToOne(2, 32);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i16.
    template<typename T>
    static Mvec<T> ei16(){
        Mvec<T> res;
        return res.componentToOne(2, 33);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1i2.
    template<typename T>
    static Mvec<T> ei1i2(){
        Mvec<T> res;
        return res.componentToOne(2, 34);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 024.
    template<typename T>
    static Mvec<T> e024(){
        Mvec<T> res;
        return res.componentToOne(2, 35);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 025.
    template<typename T>
    static Mvec<T> e025(){
        Mvec<T> res;
        return res.componentToOne(2, 36);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 026.
    template<typename T>
    static Mvec<T> e026(){
        Mvec<T> res;
        return res.componentToOne(2, 37);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 02i2.
    template<typename T>
    static Mvec<T> e02i2(){
        Mvec<T> res;
        return res.componentToOne(2, 38);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 45.
    template<typename T>
    static Mvec<T> e45(){
        Mvec<T> res;
        return res.componentToOne(2, 39);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 46.
    template<typename T>
    static Mvec<T> e46(){
        Mvec<T> res;
        return res.componentToOne(2, 40);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 4i2.
    template<typename T>
    static Mvec<T> e4i2(){
        Mvec<T> res;
        return res.componentToOne(2, 41);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 56.
    template<typename T>
    static Mvec<T> e56(){
        Mvec<T> res;
        return res.componentToOne(2, 42);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 5i2.
    template<typename T>
    static Mvec<T> e5i2(){
        Mvec<T> res;
        return res.componentToOne(2, 43);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 6i2.
    template<typename T>
    static Mvec<T> e6i2(){
        Mvec<T> res;
        return res.componentToOne(2, 44);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112.
    template<typename T>
    static Mvec<T> e0112(){
        Mvec<T> res;
        return res.componentToOne(3, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113.
    template<typename T>
    static Mvec<T> e0113(){
        Mvec<T> res;
        return res.componentToOne(3, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1.
    template<typename T>
    static Mvec<T> e011i1(){
        Mvec<T> res;
        return res.componentToOne(3, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01102.
    template<typename T>
    static Mvec<T> e01102(){
        Mvec<T> res;
        return res.componentToOne(3, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0114.
    template<typename T>
    static Mvec<T> e0114(){
        Mvec<T> res;
        return res.componentToOne(3, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0115.
    template<typename T>
    static Mvec<T> e0115(){
        Mvec<T> res;
        return res.componentToOne(3, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0116.
    template<typename T>
    static Mvec<T> e0116(){
        Mvec<T> res;
        return res.componentToOne(3, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i2.
    template<typename T>
    static Mvec<T> e011i2(){
        Mvec<T> res;
        return res.componentToOne(3, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123.
    template<typename T>
    static Mvec<T> e0123(){
        Mvec<T> res;
        return res.componentToOne(3, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1.
    template<typename T>
    static Mvec<T> e012i1(){
        Mvec<T> res;
        return res.componentToOne(3, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01202.
    template<typename T>
    static Mvec<T> e01202(){
        Mvec<T> res;
        return res.componentToOne(3, 10);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0124.
    template<typename T>
    static Mvec<T> e0124(){
        Mvec<T> res;
        return res.componentToOne(3, 11);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0125.
    template<typename T>
    static Mvec<T> e0125(){
        Mvec<T> res;
        return res.componentToOne(3, 12);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0126.
    template<typename T>
    static Mvec<T> e0126(){
        Mvec<T> res;
        return res.componentToOne(3, 13);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i2.
    template<typename T>
    static Mvec<T> e012i2(){
        Mvec<T> res;
        return res.componentToOne(3, 14);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1.
    template<typename T>
    static Mvec<T> e013i1(){
        Mvec<T> res;
        return res.componentToOne(3, 15);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01302.
    template<typename T>
    static Mvec<T> e01302(){
        Mvec<T> res;
        return res.componentToOne(3, 16);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0134.
    template<typename T>
    static Mvec<T> e0134(){
        Mvec<T> res;
        return res.componentToOne(3, 17);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0135.
    template<typename T>
    static Mvec<T> e0135(){
        Mvec<T> res;
        return res.componentToOne(3, 18);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0136.
    template<typename T>
    static Mvec<T> e0136(){
        Mvec<T> res;
        return res.componentToOne(3, 19);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i2.
    template<typename T>
    static Mvec<T> e013i2(){
        Mvec<T> res;
        return res.componentToOne(3, 20);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i102.
    template<typename T>
    static Mvec<T> e01i102(){
        Mvec<T> res;
        return res.componentToOne(3, 21);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i14.
    template<typename T>
    static Mvec<T> e01i14(){
        Mvec<T> res;
        return res.componentToOne(3, 22);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i15.
    template<typename T>
    static Mvec<T> e01i15(){
        Mvec<T> res;
        return res.componentToOne(3, 23);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i16.
    template<typename T>
    static Mvec<T> e01i16(){
        Mvec<T> res;
        return res.componentToOne(3, 24);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1i2.
    template<typename T>
    static Mvec<T> e01i1i2(){
        Mvec<T> res;
        return res.componentToOne(3, 25);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01024.
    template<typename T>
    static Mvec<T> e01024(){
        Mvec<T> res;
        return res.componentToOne(3, 26);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01025.
    template<typename T>
    static Mvec<T> e01025(){
        Mvec<T> res;
        return res.componentToOne(3, 27);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01026.
    template<typename T>
    static Mvec<T> e01026(){
        Mvec<T> res;
        return res.componentToOne(3, 28);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0102i2.
    template<typename T>
    static Mvec<T> e0102i2(){
        Mvec<T> res;
        return res.componentToOne(3, 29);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0145.
    template<typename T>
    static Mvec<T> e0145(){
        Mvec<T> res;
        return res.componentToOne(3, 30);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0146.
    template<typename T>
    static Mvec<T> e0146(){
        Mvec<T> res;
        return res.componentToOne(3, 31);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 014i2.
    template<typename T>
    static Mvec<T> e014i2(){
        Mvec<T> res;
        return res.componentToOne(3, 32);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0156.
    template<typename T>
    static Mvec<T> e0156(){
        Mvec<T> res;
        return res.componentToOne(3, 33);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 015i2.
    template<typename T>
    static Mvec<T> e015i2(){
        Mvec<T> res;
        return res.componentToOne(3, 34);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 016i2.
    template<typename T>
    static Mvec<T> e016i2(){
        Mvec<T> res;
        return res.componentToOne(3, 35);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123.
    template<typename T>
    static Mvec<T> e123(){
        Mvec<T> res;
        return res.componentToOne(3, 36);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1.
    template<typename T>
    static Mvec<T> e12i1(){
        Mvec<T> res;
        return res.componentToOne(3, 37);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1202.
    template<typename T>
    static Mvec<T> e1202(){
        Mvec<T> res;
        return res.componentToOne(3, 38);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 124.
    template<typename T>
    static Mvec<T> e124(){
        Mvec<T> res;
        return res.componentToOne(3, 39);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 125.
    template<typename T>
    static Mvec<T> e125(){
        Mvec<T> res;
        return res.componentToOne(3, 40);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 126.
    template<typename T>
    static Mvec<T> e126(){
        Mvec<T> res;
        return res.componentToOne(3, 41);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i2.
    template<typename T>
    static Mvec<T> e12i2(){
        Mvec<T> res;
        return res.componentToOne(3, 42);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1.
    template<typename T>
    static Mvec<T> e13i1(){
        Mvec<T> res;
        return res.componentToOne(3, 43);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1302.
    template<typename T>
    static Mvec<T> e1302(){
        Mvec<T> res;
        return res.componentToOne(3, 44);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 134.
    template<typename T>
    static Mvec<T> e134(){
        Mvec<T> res;
        return res.componentToOne(3, 45);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 135.
    template<typename T>
    static Mvec<T> e135(){
        Mvec<T> res;
        return res.componentToOne(3, 46);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 136.
    template<typename T>
    static Mvec<T> e136(){
        Mvec<T> res;
        return res.componentToOne(3, 47);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i2.
    template<typename T>
    static Mvec<T> e13i2(){
        Mvec<T> res;
        return res.componentToOne(3, 48);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i102.
    template<typename T>
    static Mvec<T> e1i102(){
        Mvec<T> res;
        return res.componentToOne(3, 49);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i14.
    template<typename T>
    static Mvec<T> e1i14(){
        Mvec<T> res;
        return res.componentToOne(3, 50);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i15.
    template<typename T>
    static Mvec<T> e1i15(){
        Mvec<T> res;
        return res.componentToOne(3, 51);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i16.
    template<typename T>
    static Mvec<T> e1i16(){
        Mvec<T> res;
        return res.componentToOne(3, 52);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1i2.
    template<typename T>
    static Mvec<T> e1i1i2(){
        Mvec<T> res;
        return res.componentToOne(3, 53);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1024.
    template<typename T>
    static Mvec<T> e1024(){
        Mvec<T> res;
        return res.componentToOne(3, 54);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1025.
    template<typename T>
    static Mvec<T> e1025(){
        Mvec<T> res;
        return res.componentToOne(3, 55);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1026.
    template<typename T>
    static Mvec<T> e1026(){
        Mvec<T> res;
        return res.componentToOne(3, 56);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 102i2.
    template<typename T>
    static Mvec<T> e102i2(){
        Mvec<T> res;
        return res.componentToOne(3, 57);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 145.
    template<typename T>
    static Mvec<T> e145(){
        Mvec<T> res;
        return res.componentToOne(3, 58);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 146.
    template<typename T>
    static Mvec<T> e146(){
        Mvec<T> res;
        return res.componentToOne(3, 59);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 14i2.
    template<typename T>
    static Mvec<T> e14i2(){
        Mvec<T> res;
        return res.componentToOne(3, 60);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 156.
    template<typename T>
    static Mvec<T> e156(){
        Mvec<T> res;
        return res.componentToOne(3, 61);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 15i2.
    template<typename T>
    static Mvec<T> e15i2(){
        Mvec<T> res;
        return res.componentToOne(3, 62);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 16i2.
    template<typename T>
    static Mvec<T> e16i2(){
        Mvec<T> res;
        return res.componentToOne(3, 63);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1.
    template<typename T>
    static Mvec<T> e23i1(){
        Mvec<T> res;
        return res.componentToOne(3, 64);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2302.
    template<typename T>
    static Mvec<T> e2302(){
        Mvec<T> res;
        return res.componentToOne(3, 65);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 234.
    template<typename T>
    static Mvec<T> e234(){
        Mvec<T> res;
        return res.componentToOne(3, 66);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 235.
    template<typename T>
    static Mvec<T> e235(){
        Mvec<T> res;
        return res.componentToOne(3, 67);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 236.
    template<typename T>
    static Mvec<T> e236(){
        Mvec<T> res;
        return res.componentToOne(3, 68);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i2.
    template<typename T>
    static Mvec<T> e23i2(){
        Mvec<T> res;
        return res.componentToOne(3, 69);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i102.
    template<typename T>
    static Mvec<T> e2i102(){
        Mvec<T> res;
        return res.componentToOne(3, 70);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i14.
    template<typename T>
    static Mvec<T> e2i14(){
        Mvec<T> res;
        return res.componentToOne(3, 71);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i15.
    template<typename T>
    static Mvec<T> e2i15(){
        Mvec<T> res;
        return res.componentToOne(3, 72);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i16.
    template<typename T>
    static Mvec<T> e2i16(){
        Mvec<T> res;
        return res.componentToOne(3, 73);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1i2.
    template<typename T>
    static Mvec<T> e2i1i2(){
        Mvec<T> res;
        return res.componentToOne(3, 74);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2024.
    template<typename T>
    static Mvec<T> e2024(){
        Mvec<T> res;
        return res.componentToOne(3, 75);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2025.
    template<typename T>
    static Mvec<T> e2025(){
        Mvec<T> res;
        return res.componentToOne(3, 76);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2026.
    template<typename T>
    static Mvec<T> e2026(){
        Mvec<T> res;
        return res.componentToOne(3, 77);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 202i2.
    template<typename T>
    static Mvec<T> e202i2(){
        Mvec<T> res;
        return res.componentToOne(3, 78);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 245.
    template<typename T>
    static Mvec<T> e245(){
        Mvec<T> res;
        return res.componentToOne(3, 79);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 246.
    template<typename T>
    static Mvec<T> e246(){
        Mvec<T> res;
        return res.componentToOne(3, 80);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 24i2.
    template<typename T>
    static Mvec<T> e24i2(){
        Mvec<T> res;
        return res.componentToOne(3, 81);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 256.
    template<typename T>
    static Mvec<T> e256(){
        Mvec<T> res;
        return res.componentToOne(3, 82);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 25i2.
    template<typename T>
    static Mvec<T> e25i2(){
        Mvec<T> res;
        return res.componentToOne(3, 83);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 26i2.
    template<typename T>
    static Mvec<T> e26i2(){
        Mvec<T> res;
        return res.componentToOne(3, 84);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i102.
    template<typename T>
    static Mvec<T> e3i102(){
        Mvec<T> res;
        return res.componentToOne(3, 85);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i14.
    template<typename T>
    static Mvec<T> e3i14(){
        Mvec<T> res;
        return res.componentToOne(3, 86);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i15.
    template<typename T>
    static Mvec<T> e3i15(){
        Mvec<T> res;
        return res.componentToOne(3, 87);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i16.
    template<typename T>
    static Mvec<T> e3i16(){
        Mvec<T> res;
        return res.componentToOne(3, 88);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1i2.
    template<typename T>
    static Mvec<T> e3i1i2(){
        Mvec<T> res;
        return res.componentToOne(3, 89);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3024.
    template<typename T>
    static Mvec<T> e3024(){
        Mvec<T> res;
        return res.componentToOne(3, 90);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3025.
    template<typename T>
    static Mvec<T> e3025(){
        Mvec<T> res;
        return res.componentToOne(3, 91);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3026.
    template<typename T>
    static Mvec<T> e3026(){
        Mvec<T> res;
        return res.componentToOne(3, 92);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 302i2.
    template<typename T>
    static Mvec<T> e302i2(){
        Mvec<T> res;
        return res.componentToOne(3, 93);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 345.
    template<typename T>
    static Mvec<T> e345(){
        Mvec<T> res;
        return res.componentToOne(3, 94);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 346.
    template<typename T>
    static Mvec<T> e346(){
        Mvec<T> res;
        return res.componentToOne(3, 95);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 34i2.
    template<typename T>
    static Mvec<T> e34i2(){
        Mvec<T> res;
        return res.componentToOne(3, 96);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 356.
    template<typename T>
    static Mvec<T> e356(){
        Mvec<T> res;
        return res.componentToOne(3, 97);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 35i2.
    template<typename T>
    static Mvec<T> e35i2(){
        Mvec<T> res;
        return res.componentToOne(3, 98);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 36i2.
    template<typename T>
    static Mvec<T> e36i2(){
        Mvec<T> res;
        return res.componentToOne(3, 99);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1024.
    template<typename T>
    static Mvec<T> ei1024(){
        Mvec<T> res;
        return res.componentToOne(3, 100);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1025.
    template<typename T>
    static Mvec<T> ei1025(){
        Mvec<T> res;
        return res.componentToOne(3, 101);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1026.
    template<typename T>
    static Mvec<T> ei1026(){
        Mvec<T> res;
        return res.componentToOne(3, 102);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i102i2.
    template<typename T>
    static Mvec<T> ei102i2(){
        Mvec<T> res;
        return res.componentToOne(3, 103);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i145.
    template<typename T>
    static Mvec<T> ei145(){
        Mvec<T> res;
        return res.componentToOne(3, 104);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i146.
    template<typename T>
    static Mvec<T> ei146(){
        Mvec<T> res;
        return res.componentToOne(3, 105);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i14i2.
    template<typename T>
    static Mvec<T> ei14i2(){
        Mvec<T> res;
        return res.componentToOne(3, 106);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i156.
    template<typename T>
    static Mvec<T> ei156(){
        Mvec<T> res;
        return res.componentToOne(3, 107);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i15i2.
    template<typename T>
    static Mvec<T> ei15i2(){
        Mvec<T> res;
        return res.componentToOne(3, 108);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i16i2.
    template<typename T>
    static Mvec<T> ei16i2(){
        Mvec<T> res;
        return res.componentToOne(3, 109);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0245.
    template<typename T>
    static Mvec<T> e0245(){
        Mvec<T> res;
        return res.componentToOne(3, 110);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0246.
    template<typename T>
    static Mvec<T> e0246(){
        Mvec<T> res;
        return res.componentToOne(3, 111);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 024i2.
    template<typename T>
    static Mvec<T> e024i2(){
        Mvec<T> res;
        return res.componentToOne(3, 112);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0256.
    template<typename T>
    static Mvec<T> e0256(){
        Mvec<T> res;
        return res.componentToOne(3, 113);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 025i2.
    template<typename T>
    static Mvec<T> e025i2(){
        Mvec<T> res;
        return res.componentToOne(3, 114);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 026i2.
    template<typename T>
    static Mvec<T> e026i2(){
        Mvec<T> res;
        return res.componentToOne(3, 115);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 456.
    template<typename T>
    static Mvec<T> e456(){
        Mvec<T> res;
        return res.componentToOne(3, 116);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 45i2.
    template<typename T>
    static Mvec<T> e45i2(){
        Mvec<T> res;
        return res.componentToOne(3, 117);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 46i2.
    template<typename T>
    static Mvec<T> e46i2(){
        Mvec<T> res;
        return res.componentToOne(3, 118);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 56i2.
    template<typename T>
    static Mvec<T> e56i2(){
        Mvec<T> res;
        return res.componentToOne(3, 119);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123.
    template<typename T>
    static Mvec<T> e01123(){
        Mvec<T> res;
        return res.componentToOne(4, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1.
    template<typename T>
    static Mvec<T> e0112i1(){
        Mvec<T> res;
        return res.componentToOne(4, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011202.
    template<typename T>
    static Mvec<T> e011202(){
        Mvec<T> res;
        return res.componentToOne(4, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01124.
    template<typename T>
    static Mvec<T> e01124(){
        Mvec<T> res;
        return res.componentToOne(4, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01125.
    template<typename T>
    static Mvec<T> e01125(){
        Mvec<T> res;
        return res.componentToOne(4, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01126.
    template<typename T>
    static Mvec<T> e01126(){
        Mvec<T> res;
        return res.componentToOne(4, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i2.
    template<typename T>
    static Mvec<T> e0112i2(){
        Mvec<T> res;
        return res.componentToOne(4, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1.
    template<typename T>
    static Mvec<T> e0113i1(){
        Mvec<T> res;
        return res.componentToOne(4, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011302.
    template<typename T>
    static Mvec<T> e011302(){
        Mvec<T> res;
        return res.componentToOne(4, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01134.
    template<typename T>
    static Mvec<T> e01134(){
        Mvec<T> res;
        return res.componentToOne(4, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01135.
    template<typename T>
    static Mvec<T> e01135(){
        Mvec<T> res;
        return res.componentToOne(4, 10);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01136.
    template<typename T>
    static Mvec<T> e01136(){
        Mvec<T> res;
        return res.componentToOne(4, 11);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i2.
    template<typename T>
    static Mvec<T> e0113i2(){
        Mvec<T> res;
        return res.componentToOne(4, 12);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i102.
    template<typename T>
    static Mvec<T> e011i102(){
        Mvec<T> res;
        return res.componentToOne(4, 13);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i14.
    template<typename T>
    static Mvec<T> e011i14(){
        Mvec<T> res;
        return res.componentToOne(4, 14);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i15.
    template<typename T>
    static Mvec<T> e011i15(){
        Mvec<T> res;
        return res.componentToOne(4, 15);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i16.
    template<typename T>
    static Mvec<T> e011i16(){
        Mvec<T> res;
        return res.componentToOne(4, 16);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1i2.
    template<typename T>
    static Mvec<T> e011i1i2(){
        Mvec<T> res;
        return res.componentToOne(4, 17);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011024.
    template<typename T>
    static Mvec<T> e011024(){
        Mvec<T> res;
        return res.componentToOne(4, 18);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011025.
    template<typename T>
    static Mvec<T> e011025(){
        Mvec<T> res;
        return res.componentToOne(4, 19);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011026.
    template<typename T>
    static Mvec<T> e011026(){
        Mvec<T> res;
        return res.componentToOne(4, 20);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01102i2.
    template<typename T>
    static Mvec<T> e01102i2(){
        Mvec<T> res;
        return res.componentToOne(4, 21);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01145.
    template<typename T>
    static Mvec<T> e01145(){
        Mvec<T> res;
        return res.componentToOne(4, 22);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01146.
    template<typename T>
    static Mvec<T> e01146(){
        Mvec<T> res;
        return res.componentToOne(4, 23);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0114i2.
    template<typename T>
    static Mvec<T> e0114i2(){
        Mvec<T> res;
        return res.componentToOne(4, 24);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01156.
    template<typename T>
    static Mvec<T> e01156(){
        Mvec<T> res;
        return res.componentToOne(4, 25);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0115i2.
    template<typename T>
    static Mvec<T> e0115i2(){
        Mvec<T> res;
        return res.componentToOne(4, 26);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0116i2.
    template<typename T>
    static Mvec<T> e0116i2(){
        Mvec<T> res;
        return res.componentToOne(4, 27);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1.
    template<typename T>
    static Mvec<T> e0123i1(){
        Mvec<T> res;
        return res.componentToOne(4, 28);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012302.
    template<typename T>
    static Mvec<T> e012302(){
        Mvec<T> res;
        return res.componentToOne(4, 29);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01234.
    template<typename T>
    static Mvec<T> e01234(){
        Mvec<T> res;
        return res.componentToOne(4, 30);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01235.
    template<typename T>
    static Mvec<T> e01235(){
        Mvec<T> res;
        return res.componentToOne(4, 31);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01236.
    template<typename T>
    static Mvec<T> e01236(){
        Mvec<T> res;
        return res.componentToOne(4, 32);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i2.
    template<typename T>
    static Mvec<T> e0123i2(){
        Mvec<T> res;
        return res.componentToOne(4, 33);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i102.
    template<typename T>
    static Mvec<T> e012i102(){
        Mvec<T> res;
        return res.componentToOne(4, 34);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i14.
    template<typename T>
    static Mvec<T> e012i14(){
        Mvec<T> res;
        return res.componentToOne(4, 35);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i15.
    template<typename T>
    static Mvec<T> e012i15(){
        Mvec<T> res;
        return res.componentToOne(4, 36);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i16.
    template<typename T>
    static Mvec<T> e012i16(){
        Mvec<T> res;
        return res.componentToOne(4, 37);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1i2.
    template<typename T>
    static Mvec<T> e012i1i2(){
        Mvec<T> res;
        return res.componentToOne(4, 38);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012024.
    template<typename T>
    static Mvec<T> e012024(){
        Mvec<T> res;
        return res.componentToOne(4, 39);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012025.
    template<typename T>
    static Mvec<T> e012025(){
        Mvec<T> res;
        return res.componentToOne(4, 40);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012026.
    template<typename T>
    static Mvec<T> e012026(){
        Mvec<T> res;
        return res.componentToOne(4, 41);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01202i2.
    template<typename T>
    static Mvec<T> e01202i2(){
        Mvec<T> res;
        return res.componentToOne(4, 42);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01245.
    template<typename T>
    static Mvec<T> e01245(){
        Mvec<T> res;
        return res.componentToOne(4, 43);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01246.
    template<typename T>
    static Mvec<T> e01246(){
        Mvec<T> res;
        return res.componentToOne(4, 44);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0124i2.
    template<typename T>
    static Mvec<T> e0124i2(){
        Mvec<T> res;
        return res.componentToOne(4, 45);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01256.
    template<typename T>
    static Mvec<T> e01256(){
        Mvec<T> res;
        return res.componentToOne(4, 46);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0125i2.
    template<typename T>
    static Mvec<T> e0125i2(){
        Mvec<T> res;
        return res.componentToOne(4, 47);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0126i2.
    template<typename T>
    static Mvec<T> e0126i2(){
        Mvec<T> res;
        return res.componentToOne(4, 48);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i102.
    template<typename T>
    static Mvec<T> e013i102(){
        Mvec<T> res;
        return res.componentToOne(4, 49);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i14.
    template<typename T>
    static Mvec<T> e013i14(){
        Mvec<T> res;
        return res.componentToOne(4, 50);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i15.
    template<typename T>
    static Mvec<T> e013i15(){
        Mvec<T> res;
        return res.componentToOne(4, 51);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i16.
    template<typename T>
    static Mvec<T> e013i16(){
        Mvec<T> res;
        return res.componentToOne(4, 52);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1i2.
    template<typename T>
    static Mvec<T> e013i1i2(){
        Mvec<T> res;
        return res.componentToOne(4, 53);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013024.
    template<typename T>
    static Mvec<T> e013024(){
        Mvec<T> res;
        return res.componentToOne(4, 54);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013025.
    template<typename T>
    static Mvec<T> e013025(){
        Mvec<T> res;
        return res.componentToOne(4, 55);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013026.
    template<typename T>
    static Mvec<T> e013026(){
        Mvec<T> res;
        return res.componentToOne(4, 56);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01302i2.
    template<typename T>
    static Mvec<T> e01302i2(){
        Mvec<T> res;
        return res.componentToOne(4, 57);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01345.
    template<typename T>
    static Mvec<T> e01345(){
        Mvec<T> res;
        return res.componentToOne(4, 58);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01346.
    template<typename T>
    static Mvec<T> e01346(){
        Mvec<T> res;
        return res.componentToOne(4, 59);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0134i2.
    template<typename T>
    static Mvec<T> e0134i2(){
        Mvec<T> res;
        return res.componentToOne(4, 60);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01356.
    template<typename T>
    static Mvec<T> e01356(){
        Mvec<T> res;
        return res.componentToOne(4, 61);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0135i2.
    template<typename T>
    static Mvec<T> e0135i2(){
        Mvec<T> res;
        return res.componentToOne(4, 62);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0136i2.
    template<typename T>
    static Mvec<T> e0136i2(){
        Mvec<T> res;
        return res.componentToOne(4, 63);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1024.
    template<typename T>
    static Mvec<T> e01i1024(){
        Mvec<T> res;
        return res.componentToOne(4, 64);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1025.
    template<typename T>
    static Mvec<T> e01i1025(){
        Mvec<T> res;
        return res.componentToOne(4, 65);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1026.
    template<typename T>
    static Mvec<T> e01i1026(){
        Mvec<T> res;
        return res.componentToOne(4, 66);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i102i2.
    template<typename T>
    static Mvec<T> e01i102i2(){
        Mvec<T> res;
        return res.componentToOne(4, 67);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i145.
    template<typename T>
    static Mvec<T> e01i145(){
        Mvec<T> res;
        return res.componentToOne(4, 68);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i146.
    template<typename T>
    static Mvec<T> e01i146(){
        Mvec<T> res;
        return res.componentToOne(4, 69);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i14i2.
    template<typename T>
    static Mvec<T> e01i14i2(){
        Mvec<T> res;
        return res.componentToOne(4, 70);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i156.
    template<typename T>
    static Mvec<T> e01i156(){
        Mvec<T> res;
        return res.componentToOne(4, 71);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i15i2.
    template<typename T>
    static Mvec<T> e01i15i2(){
        Mvec<T> res;
        return res.componentToOne(4, 72);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i16i2.
    template<typename T>
    static Mvec<T> e01i16i2(){
        Mvec<T> res;
        return res.componentToOne(4, 73);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 010245.
    template<typename T>
    static Mvec<T> e010245(){
        Mvec<T> res;
        return res.componentToOne(4, 74);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 010246.
    template<typename T>
    static Mvec<T> e010246(){
        Mvec<T> res;
        return res.componentToOne(4, 75);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01024i2.
    template<typename T>
    static Mvec<T> e01024i2(){
        Mvec<T> res;
        return res.componentToOne(4, 76);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 010256.
    template<typename T>
    static Mvec<T> e010256(){
        Mvec<T> res;
        return res.componentToOne(4, 77);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01025i2.
    template<typename T>
    static Mvec<T> e01025i2(){
        Mvec<T> res;
        return res.componentToOne(4, 78);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01026i2.
    template<typename T>
    static Mvec<T> e01026i2(){
        Mvec<T> res;
        return res.componentToOne(4, 79);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01456.
    template<typename T>
    static Mvec<T> e01456(){
        Mvec<T> res;
        return res.componentToOne(4, 80);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0145i2.
    template<typename T>
    static Mvec<T> e0145i2(){
        Mvec<T> res;
        return res.componentToOne(4, 81);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0146i2.
    template<typename T>
    static Mvec<T> e0146i2(){
        Mvec<T> res;
        return res.componentToOne(4, 82);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0156i2.
    template<typename T>
    static Mvec<T> e0156i2(){
        Mvec<T> res;
        return res.componentToOne(4, 83);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1.
    template<typename T>
    static Mvec<T> e123i1(){
        Mvec<T> res;
        return res.componentToOne(4, 84);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12302.
    template<typename T>
    static Mvec<T> e12302(){
        Mvec<T> res;
        return res.componentToOne(4, 85);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1234.
    template<typename T>
    static Mvec<T> e1234(){
        Mvec<T> res;
        return res.componentToOne(4, 86);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1235.
    template<typename T>
    static Mvec<T> e1235(){
        Mvec<T> res;
        return res.componentToOne(4, 87);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1236.
    template<typename T>
    static Mvec<T> e1236(){
        Mvec<T> res;
        return res.componentToOne(4, 88);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i2.
    template<typename T>
    static Mvec<T> e123i2(){
        Mvec<T> res;
        return res.componentToOne(4, 89);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i102.
    template<typename T>
    static Mvec<T> e12i102(){
        Mvec<T> res;
        return res.componentToOne(4, 90);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i14.
    template<typename T>
    static Mvec<T> e12i14(){
        Mvec<T> res;
        return res.componentToOne(4, 91);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i15.
    template<typename T>
    static Mvec<T> e12i15(){
        Mvec<T> res;
        return res.componentToOne(4, 92);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i16.
    template<typename T>
    static Mvec<T> e12i16(){
        Mvec<T> res;
        return res.componentToOne(4, 93);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1i2.
    template<typename T>
    static Mvec<T> e12i1i2(){
        Mvec<T> res;
        return res.componentToOne(4, 94);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12024.
    template<typename T>
    static Mvec<T> e12024(){
        Mvec<T> res;
        return res.componentToOne(4, 95);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12025.
    template<typename T>
    static Mvec<T> e12025(){
        Mvec<T> res;
        return res.componentToOne(4, 96);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12026.
    template<typename T>
    static Mvec<T> e12026(){
        Mvec<T> res;
        return res.componentToOne(4, 97);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1202i2.
    template<typename T>
    static Mvec<T> e1202i2(){
        Mvec<T> res;
        return res.componentToOne(4, 98);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1245.
    template<typename T>
    static Mvec<T> e1245(){
        Mvec<T> res;
        return res.componentToOne(4, 99);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1246.
    template<typename T>
    static Mvec<T> e1246(){
        Mvec<T> res;
        return res.componentToOne(4, 100);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 124i2.
    template<typename T>
    static Mvec<T> e124i2(){
        Mvec<T> res;
        return res.componentToOne(4, 101);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1256.
    template<typename T>
    static Mvec<T> e1256(){
        Mvec<T> res;
        return res.componentToOne(4, 102);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 125i2.
    template<typename T>
    static Mvec<T> e125i2(){
        Mvec<T> res;
        return res.componentToOne(4, 103);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 126i2.
    template<typename T>
    static Mvec<T> e126i2(){
        Mvec<T> res;
        return res.componentToOne(4, 104);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i102.
    template<typename T>
    static Mvec<T> e13i102(){
        Mvec<T> res;
        return res.componentToOne(4, 105);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i14.
    template<typename T>
    static Mvec<T> e13i14(){
        Mvec<T> res;
        return res.componentToOne(4, 106);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i15.
    template<typename T>
    static Mvec<T> e13i15(){
        Mvec<T> res;
        return res.componentToOne(4, 107);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i16.
    template<typename T>
    static Mvec<T> e13i16(){
        Mvec<T> res;
        return res.componentToOne(4, 108);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1i2.
    template<typename T>
    static Mvec<T> e13i1i2(){
        Mvec<T> res;
        return res.componentToOne(4, 109);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13024.
    template<typename T>
    static Mvec<T> e13024(){
        Mvec<T> res;
        return res.componentToOne(4, 110);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13025.
    template<typename T>
    static Mvec<T> e13025(){
        Mvec<T> res;
        return res.componentToOne(4, 111);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13026.
    template<typename T>
    static Mvec<T> e13026(){
        Mvec<T> res;
        return res.componentToOne(4, 112);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1302i2.
    template<typename T>
    static Mvec<T> e1302i2(){
        Mvec<T> res;
        return res.componentToOne(4, 113);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1345.
    template<typename T>
    static Mvec<T> e1345(){
        Mvec<T> res;
        return res.componentToOne(4, 114);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1346.
    template<typename T>
    static Mvec<T> e1346(){
        Mvec<T> res;
        return res.componentToOne(4, 115);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 134i2.
    template<typename T>
    static Mvec<T> e134i2(){
        Mvec<T> res;
        return res.componentToOne(4, 116);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1356.
    template<typename T>
    static Mvec<T> e1356(){
        Mvec<T> res;
        return res.componentToOne(4, 117);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 135i2.
    template<typename T>
    static Mvec<T> e135i2(){
        Mvec<T> res;
        return res.componentToOne(4, 118);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 136i2.
    template<typename T>
    static Mvec<T> e136i2(){
        Mvec<T> res;
        return res.componentToOne(4, 119);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1024.
    template<typename T>
    static Mvec<T> e1i1024(){
        Mvec<T> res;
        return res.componentToOne(4, 120);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1025.
    template<typename T>
    static Mvec<T> e1i1025(){
        Mvec<T> res;
        return res.componentToOne(4, 121);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1026.
    template<typename T>
    static Mvec<T> e1i1026(){
        Mvec<T> res;
        return res.componentToOne(4, 122);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i102i2.
    template<typename T>
    static Mvec<T> e1i102i2(){
        Mvec<T> res;
        return res.componentToOne(4, 123);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i145.
    template<typename T>
    static Mvec<T> e1i145(){
        Mvec<T> res;
        return res.componentToOne(4, 124);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i146.
    template<typename T>
    static Mvec<T> e1i146(){
        Mvec<T> res;
        return res.componentToOne(4, 125);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i14i2.
    template<typename T>
    static Mvec<T> e1i14i2(){
        Mvec<T> res;
        return res.componentToOne(4, 126);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i156.
    template<typename T>
    static Mvec<T> e1i156(){
        Mvec<T> res;
        return res.componentToOne(4, 127);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i15i2.
    template<typename T>
    static Mvec<T> e1i15i2(){
        Mvec<T> res;
        return res.componentToOne(4, 128);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i16i2.
    template<typename T>
    static Mvec<T> e1i16i2(){
        Mvec<T> res;
        return res.componentToOne(4, 129);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 10245.
    template<typename T>
    static Mvec<T> e10245(){
        Mvec<T> res;
        return res.componentToOne(4, 130);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 10246.
    template<typename T>
    static Mvec<T> e10246(){
        Mvec<T> res;
        return res.componentToOne(4, 131);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1024i2.
    template<typename T>
    static Mvec<T> e1024i2(){
        Mvec<T> res;
        return res.componentToOne(4, 132);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 10256.
    template<typename T>
    static Mvec<T> e10256(){
        Mvec<T> res;
        return res.componentToOne(4, 133);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1025i2.
    template<typename T>
    static Mvec<T> e1025i2(){
        Mvec<T> res;
        return res.componentToOne(4, 134);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1026i2.
    template<typename T>
    static Mvec<T> e1026i2(){
        Mvec<T> res;
        return res.componentToOne(4, 135);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1456.
    template<typename T>
    static Mvec<T> e1456(){
        Mvec<T> res;
        return res.componentToOne(4, 136);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 145i2.
    template<typename T>
    static Mvec<T> e145i2(){
        Mvec<T> res;
        return res.componentToOne(4, 137);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 146i2.
    template<typename T>
    static Mvec<T> e146i2(){
        Mvec<T> res;
        return res.componentToOne(4, 138);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 156i2.
    template<typename T>
    static Mvec<T> e156i2(){
        Mvec<T> res;
        return res.componentToOne(4, 139);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i102.
    template<typename T>
    static Mvec<T> e23i102(){
        Mvec<T> res;
        return res.componentToOne(4, 140);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i14.
    template<typename T>
    static Mvec<T> e23i14(){
        Mvec<T> res;
        return res.componentToOne(4, 141);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i15.
    template<typename T>
    static Mvec<T> e23i15(){
        Mvec<T> res;
        return res.componentToOne(4, 142);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i16.
    template<typename T>
    static Mvec<T> e23i16(){
        Mvec<T> res;
        return res.componentToOne(4, 143);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1i2.
    template<typename T>
    static Mvec<T> e23i1i2(){
        Mvec<T> res;
        return res.componentToOne(4, 144);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23024.
    template<typename T>
    static Mvec<T> e23024(){
        Mvec<T> res;
        return res.componentToOne(4, 145);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23025.
    template<typename T>
    static Mvec<T> e23025(){
        Mvec<T> res;
        return res.componentToOne(4, 146);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23026.
    template<typename T>
    static Mvec<T> e23026(){
        Mvec<T> res;
        return res.componentToOne(4, 147);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2302i2.
    template<typename T>
    static Mvec<T> e2302i2(){
        Mvec<T> res;
        return res.componentToOne(4, 148);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2345.
    template<typename T>
    static Mvec<T> e2345(){
        Mvec<T> res;
        return res.componentToOne(4, 149);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2346.
    template<typename T>
    static Mvec<T> e2346(){
        Mvec<T> res;
        return res.componentToOne(4, 150);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 234i2.
    template<typename T>
    static Mvec<T> e234i2(){
        Mvec<T> res;
        return res.componentToOne(4, 151);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2356.
    template<typename T>
    static Mvec<T> e2356(){
        Mvec<T> res;
        return res.componentToOne(4, 152);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 235i2.
    template<typename T>
    static Mvec<T> e235i2(){
        Mvec<T> res;
        return res.componentToOne(4, 153);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 236i2.
    template<typename T>
    static Mvec<T> e236i2(){
        Mvec<T> res;
        return res.componentToOne(4, 154);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1024.
    template<typename T>
    static Mvec<T> e2i1024(){
        Mvec<T> res;
        return res.componentToOne(4, 155);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1025.
    template<typename T>
    static Mvec<T> e2i1025(){
        Mvec<T> res;
        return res.componentToOne(4, 156);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1026.
    template<typename T>
    static Mvec<T> e2i1026(){
        Mvec<T> res;
        return res.componentToOne(4, 157);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i102i2.
    template<typename T>
    static Mvec<T> e2i102i2(){
        Mvec<T> res;
        return res.componentToOne(4, 158);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i145.
    template<typename T>
    static Mvec<T> e2i145(){
        Mvec<T> res;
        return res.componentToOne(4, 159);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i146.
    template<typename T>
    static Mvec<T> e2i146(){
        Mvec<T> res;
        return res.componentToOne(4, 160);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i14i2.
    template<typename T>
    static Mvec<T> e2i14i2(){
        Mvec<T> res;
        return res.componentToOne(4, 161);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i156.
    template<typename T>
    static Mvec<T> e2i156(){
        Mvec<T> res;
        return res.componentToOne(4, 162);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i15i2.
    template<typename T>
    static Mvec<T> e2i15i2(){
        Mvec<T> res;
        return res.componentToOne(4, 163);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i16i2.
    template<typename T>
    static Mvec<T> e2i16i2(){
        Mvec<T> res;
        return res.componentToOne(4, 164);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 20245.
    template<typename T>
    static Mvec<T> e20245(){
        Mvec<T> res;
        return res.componentToOne(4, 165);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 20246.
    template<typename T>
    static Mvec<T> e20246(){
        Mvec<T> res;
        return res.componentToOne(4, 166);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2024i2.
    template<typename T>
    static Mvec<T> e2024i2(){
        Mvec<T> res;
        return res.componentToOne(4, 167);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 20256.
    template<typename T>
    static Mvec<T> e20256(){
        Mvec<T> res;
        return res.componentToOne(4, 168);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2025i2.
    template<typename T>
    static Mvec<T> e2025i2(){
        Mvec<T> res;
        return res.componentToOne(4, 169);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2026i2.
    template<typename T>
    static Mvec<T> e2026i2(){
        Mvec<T> res;
        return res.componentToOne(4, 170);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2456.
    template<typename T>
    static Mvec<T> e2456(){
        Mvec<T> res;
        return res.componentToOne(4, 171);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 245i2.
    template<typename T>
    static Mvec<T> e245i2(){
        Mvec<T> res;
        return res.componentToOne(4, 172);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 246i2.
    template<typename T>
    static Mvec<T> e246i2(){
        Mvec<T> res;
        return res.componentToOne(4, 173);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 256i2.
    template<typename T>
    static Mvec<T> e256i2(){
        Mvec<T> res;
        return res.componentToOne(4, 174);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1024.
    template<typename T>
    static Mvec<T> e3i1024(){
        Mvec<T> res;
        return res.componentToOne(4, 175);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1025.
    template<typename T>
    static Mvec<T> e3i1025(){
        Mvec<T> res;
        return res.componentToOne(4, 176);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1026.
    template<typename T>
    static Mvec<T> e3i1026(){
        Mvec<T> res;
        return res.componentToOne(4, 177);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i102i2.
    template<typename T>
    static Mvec<T> e3i102i2(){
        Mvec<T> res;
        return res.componentToOne(4, 178);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i145.
    template<typename T>
    static Mvec<T> e3i145(){
        Mvec<T> res;
        return res.componentToOne(4, 179);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i146.
    template<typename T>
    static Mvec<T> e3i146(){
        Mvec<T> res;
        return res.componentToOne(4, 180);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i14i2.
    template<typename T>
    static Mvec<T> e3i14i2(){
        Mvec<T> res;
        return res.componentToOne(4, 181);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i156.
    template<typename T>
    static Mvec<T> e3i156(){
        Mvec<T> res;
        return res.componentToOne(4, 182);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i15i2.
    template<typename T>
    static Mvec<T> e3i15i2(){
        Mvec<T> res;
        return res.componentToOne(4, 183);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i16i2.
    template<typename T>
    static Mvec<T> e3i16i2(){
        Mvec<T> res;
        return res.componentToOne(4, 184);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 30245.
    template<typename T>
    static Mvec<T> e30245(){
        Mvec<T> res;
        return res.componentToOne(4, 185);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 30246.
    template<typename T>
    static Mvec<T> e30246(){
        Mvec<T> res;
        return res.componentToOne(4, 186);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3024i2.
    template<typename T>
    static Mvec<T> e3024i2(){
        Mvec<T> res;
        return res.componentToOne(4, 187);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 30256.
    template<typename T>
    static Mvec<T> e30256(){
        Mvec<T> res;
        return res.componentToOne(4, 188);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3025i2.
    template<typename T>
    static Mvec<T> e3025i2(){
        Mvec<T> res;
        return res.componentToOne(4, 189);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3026i2.
    template<typename T>
    static Mvec<T> e3026i2(){
        Mvec<T> res;
        return res.componentToOne(4, 190);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3456.
    template<typename T>
    static Mvec<T> e3456(){
        Mvec<T> res;
        return res.componentToOne(4, 191);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 345i2.
    template<typename T>
    static Mvec<T> e345i2(){
        Mvec<T> res;
        return res.componentToOne(4, 192);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 346i2.
    template<typename T>
    static Mvec<T> e346i2(){
        Mvec<T> res;
        return res.componentToOne(4, 193);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 356i2.
    template<typename T>
    static Mvec<T> e356i2(){
        Mvec<T> res;
        return res.componentToOne(4, 194);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i10245.
    template<typename T>
    static Mvec<T> ei10245(){
        Mvec<T> res;
        return res.componentToOne(4, 195);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i10246.
    template<typename T>
    static Mvec<T> ei10246(){
        Mvec<T> res;
        return res.componentToOne(4, 196);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1024i2.
    template<typename T>
    static Mvec<T> ei1024i2(){
        Mvec<T> res;
        return res.componentToOne(4, 197);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i10256.
    template<typename T>
    static Mvec<T> ei10256(){
        Mvec<T> res;
        return res.componentToOne(4, 198);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1025i2.
    template<typename T>
    static Mvec<T> ei1025i2(){
        Mvec<T> res;
        return res.componentToOne(4, 199);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1026i2.
    template<typename T>
    static Mvec<T> ei1026i2(){
        Mvec<T> res;
        return res.componentToOne(4, 200);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1456.
    template<typename T>
    static Mvec<T> ei1456(){
        Mvec<T> res;
        return res.componentToOne(4, 201);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i145i2.
    template<typename T>
    static Mvec<T> ei145i2(){
        Mvec<T> res;
        return res.componentToOne(4, 202);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i146i2.
    template<typename T>
    static Mvec<T> ei146i2(){
        Mvec<T> res;
        return res.componentToOne(4, 203);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i156i2.
    template<typename T>
    static Mvec<T> ei156i2(){
        Mvec<T> res;
        return res.componentToOne(4, 204);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 02456.
    template<typename T>
    static Mvec<T> e02456(){
        Mvec<T> res;
        return res.componentToOne(4, 205);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0245i2.
    template<typename T>
    static Mvec<T> e0245i2(){
        Mvec<T> res;
        return res.componentToOne(4, 206);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0246i2.
    template<typename T>
    static Mvec<T> e0246i2(){
        Mvec<T> res;
        return res.componentToOne(4, 207);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0256i2.
    template<typename T>
    static Mvec<T> e0256i2(){
        Mvec<T> res;
        return res.componentToOne(4, 208);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 456i2.
    template<typename T>
    static Mvec<T> e456i2(){
        Mvec<T> res;
        return res.componentToOne(4, 209);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1.
    template<typename T>
    static Mvec<T> e01123i1(){
        Mvec<T> res;
        return res.componentToOne(5, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112302.
    template<typename T>
    static Mvec<T> e0112302(){
        Mvec<T> res;
        return res.componentToOne(5, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011234.
    template<typename T>
    static Mvec<T> e011234(){
        Mvec<T> res;
        return res.componentToOne(5, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011235.
    template<typename T>
    static Mvec<T> e011235(){
        Mvec<T> res;
        return res.componentToOne(5, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011236.
    template<typename T>
    static Mvec<T> e011236(){
        Mvec<T> res;
        return res.componentToOne(5, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i2.
    template<typename T>
    static Mvec<T> e01123i2(){
        Mvec<T> res;
        return res.componentToOne(5, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i102.
    template<typename T>
    static Mvec<T> e0112i102(){
        Mvec<T> res;
        return res.componentToOne(5, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i14.
    template<typename T>
    static Mvec<T> e0112i14(){
        Mvec<T> res;
        return res.componentToOne(5, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i15.
    template<typename T>
    static Mvec<T> e0112i15(){
        Mvec<T> res;
        return res.componentToOne(5, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i16.
    template<typename T>
    static Mvec<T> e0112i16(){
        Mvec<T> res;
        return res.componentToOne(5, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1i2.
    template<typename T>
    static Mvec<T> e0112i1i2(){
        Mvec<T> res;
        return res.componentToOne(5, 10);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112024.
    template<typename T>
    static Mvec<T> e0112024(){
        Mvec<T> res;
        return res.componentToOne(5, 11);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112025.
    template<typename T>
    static Mvec<T> e0112025(){
        Mvec<T> res;
        return res.componentToOne(5, 12);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112026.
    template<typename T>
    static Mvec<T> e0112026(){
        Mvec<T> res;
        return res.componentToOne(5, 13);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011202i2.
    template<typename T>
    static Mvec<T> e011202i2(){
        Mvec<T> res;
        return res.componentToOne(5, 14);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011245.
    template<typename T>
    static Mvec<T> e011245(){
        Mvec<T> res;
        return res.componentToOne(5, 15);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011246.
    template<typename T>
    static Mvec<T> e011246(){
        Mvec<T> res;
        return res.componentToOne(5, 16);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01124i2.
    template<typename T>
    static Mvec<T> e01124i2(){
        Mvec<T> res;
        return res.componentToOne(5, 17);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011256.
    template<typename T>
    static Mvec<T> e011256(){
        Mvec<T> res;
        return res.componentToOne(5, 18);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01125i2.
    template<typename T>
    static Mvec<T> e01125i2(){
        Mvec<T> res;
        return res.componentToOne(5, 19);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01126i2.
    template<typename T>
    static Mvec<T> e01126i2(){
        Mvec<T> res;
        return res.componentToOne(5, 20);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i102.
    template<typename T>
    static Mvec<T> e0113i102(){
        Mvec<T> res;
        return res.componentToOne(5, 21);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i14.
    template<typename T>
    static Mvec<T> e0113i14(){
        Mvec<T> res;
        return res.componentToOne(5, 22);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i15.
    template<typename T>
    static Mvec<T> e0113i15(){
        Mvec<T> res;
        return res.componentToOne(5, 23);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i16.
    template<typename T>
    static Mvec<T> e0113i16(){
        Mvec<T> res;
        return res.componentToOne(5, 24);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1i2.
    template<typename T>
    static Mvec<T> e0113i1i2(){
        Mvec<T> res;
        return res.componentToOne(5, 25);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113024.
    template<typename T>
    static Mvec<T> e0113024(){
        Mvec<T> res;
        return res.componentToOne(5, 26);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113025.
    template<typename T>
    static Mvec<T> e0113025(){
        Mvec<T> res;
        return res.componentToOne(5, 27);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113026.
    template<typename T>
    static Mvec<T> e0113026(){
        Mvec<T> res;
        return res.componentToOne(5, 28);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011302i2.
    template<typename T>
    static Mvec<T> e011302i2(){
        Mvec<T> res;
        return res.componentToOne(5, 29);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011345.
    template<typename T>
    static Mvec<T> e011345(){
        Mvec<T> res;
        return res.componentToOne(5, 30);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011346.
    template<typename T>
    static Mvec<T> e011346(){
        Mvec<T> res;
        return res.componentToOne(5, 31);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01134i2.
    template<typename T>
    static Mvec<T> e01134i2(){
        Mvec<T> res;
        return res.componentToOne(5, 32);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011356.
    template<typename T>
    static Mvec<T> e011356(){
        Mvec<T> res;
        return res.componentToOne(5, 33);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01135i2.
    template<typename T>
    static Mvec<T> e01135i2(){
        Mvec<T> res;
        return res.componentToOne(5, 34);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01136i2.
    template<typename T>
    static Mvec<T> e01136i2(){
        Mvec<T> res;
        return res.componentToOne(5, 35);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1024.
    template<typename T>
    static Mvec<T> e011i1024(){
        Mvec<T> res;
        return res.componentToOne(5, 36);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1025.
    template<typename T>
    static Mvec<T> e011i1025(){
        Mvec<T> res;
        return res.componentToOne(5, 37);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1026.
    template<typename T>
    static Mvec<T> e011i1026(){
        Mvec<T> res;
        return res.componentToOne(5, 38);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i102i2.
    template<typename T>
    static Mvec<T> e011i102i2(){
        Mvec<T> res;
        return res.componentToOne(5, 39);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i145.
    template<typename T>
    static Mvec<T> e011i145(){
        Mvec<T> res;
        return res.componentToOne(5, 40);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i146.
    template<typename T>
    static Mvec<T> e011i146(){
        Mvec<T> res;
        return res.componentToOne(5, 41);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i14i2.
    template<typename T>
    static Mvec<T> e011i14i2(){
        Mvec<T> res;
        return res.componentToOne(5, 42);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i156.
    template<typename T>
    static Mvec<T> e011i156(){
        Mvec<T> res;
        return res.componentToOne(5, 43);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i15i2.
    template<typename T>
    static Mvec<T> e011i15i2(){
        Mvec<T> res;
        return res.componentToOne(5, 44);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i16i2.
    template<typename T>
    static Mvec<T> e011i16i2(){
        Mvec<T> res;
        return res.componentToOne(5, 45);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0110245.
    template<typename T>
    static Mvec<T> e0110245(){
        Mvec<T> res;
        return res.componentToOne(5, 46);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0110246.
    template<typename T>
    static Mvec<T> e0110246(){
        Mvec<T> res;
        return res.componentToOne(5, 47);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011024i2.
    template<typename T>
    static Mvec<T> e011024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 48);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0110256.
    template<typename T>
    static Mvec<T> e0110256(){
        Mvec<T> res;
        return res.componentToOne(5, 49);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011025i2.
    template<typename T>
    static Mvec<T> e011025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 50);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011026i2.
    template<typename T>
    static Mvec<T> e011026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 51);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011456.
    template<typename T>
    static Mvec<T> e011456(){
        Mvec<T> res;
        return res.componentToOne(5, 52);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01145i2.
    template<typename T>
    static Mvec<T> e01145i2(){
        Mvec<T> res;
        return res.componentToOne(5, 53);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01146i2.
    template<typename T>
    static Mvec<T> e01146i2(){
        Mvec<T> res;
        return res.componentToOne(5, 54);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01156i2.
    template<typename T>
    static Mvec<T> e01156i2(){
        Mvec<T> res;
        return res.componentToOne(5, 55);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i102.
    template<typename T>
    static Mvec<T> e0123i102(){
        Mvec<T> res;
        return res.componentToOne(5, 56);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i14.
    template<typename T>
    static Mvec<T> e0123i14(){
        Mvec<T> res;
        return res.componentToOne(5, 57);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i15.
    template<typename T>
    static Mvec<T> e0123i15(){
        Mvec<T> res;
        return res.componentToOne(5, 58);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i16.
    template<typename T>
    static Mvec<T> e0123i16(){
        Mvec<T> res;
        return res.componentToOne(5, 59);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1i2.
    template<typename T>
    static Mvec<T> e0123i1i2(){
        Mvec<T> res;
        return res.componentToOne(5, 60);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123024.
    template<typename T>
    static Mvec<T> e0123024(){
        Mvec<T> res;
        return res.componentToOne(5, 61);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123025.
    template<typename T>
    static Mvec<T> e0123025(){
        Mvec<T> res;
        return res.componentToOne(5, 62);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123026.
    template<typename T>
    static Mvec<T> e0123026(){
        Mvec<T> res;
        return res.componentToOne(5, 63);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012302i2.
    template<typename T>
    static Mvec<T> e012302i2(){
        Mvec<T> res;
        return res.componentToOne(5, 64);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012345.
    template<typename T>
    static Mvec<T> e012345(){
        Mvec<T> res;
        return res.componentToOne(5, 65);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012346.
    template<typename T>
    static Mvec<T> e012346(){
        Mvec<T> res;
        return res.componentToOne(5, 66);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01234i2.
    template<typename T>
    static Mvec<T> e01234i2(){
        Mvec<T> res;
        return res.componentToOne(5, 67);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012356.
    template<typename T>
    static Mvec<T> e012356(){
        Mvec<T> res;
        return res.componentToOne(5, 68);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01235i2.
    template<typename T>
    static Mvec<T> e01235i2(){
        Mvec<T> res;
        return res.componentToOne(5, 69);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01236i2.
    template<typename T>
    static Mvec<T> e01236i2(){
        Mvec<T> res;
        return res.componentToOne(5, 70);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1024.
    template<typename T>
    static Mvec<T> e012i1024(){
        Mvec<T> res;
        return res.componentToOne(5, 71);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1025.
    template<typename T>
    static Mvec<T> e012i1025(){
        Mvec<T> res;
        return res.componentToOne(5, 72);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1026.
    template<typename T>
    static Mvec<T> e012i1026(){
        Mvec<T> res;
        return res.componentToOne(5, 73);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i102i2.
    template<typename T>
    static Mvec<T> e012i102i2(){
        Mvec<T> res;
        return res.componentToOne(5, 74);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i145.
    template<typename T>
    static Mvec<T> e012i145(){
        Mvec<T> res;
        return res.componentToOne(5, 75);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i146.
    template<typename T>
    static Mvec<T> e012i146(){
        Mvec<T> res;
        return res.componentToOne(5, 76);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i14i2.
    template<typename T>
    static Mvec<T> e012i14i2(){
        Mvec<T> res;
        return res.componentToOne(5, 77);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i156.
    template<typename T>
    static Mvec<T> e012i156(){
        Mvec<T> res;
        return res.componentToOne(5, 78);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i15i2.
    template<typename T>
    static Mvec<T> e012i15i2(){
        Mvec<T> res;
        return res.componentToOne(5, 79);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i16i2.
    template<typename T>
    static Mvec<T> e012i16i2(){
        Mvec<T> res;
        return res.componentToOne(5, 80);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0120245.
    template<typename T>
    static Mvec<T> e0120245(){
        Mvec<T> res;
        return res.componentToOne(5, 81);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0120246.
    template<typename T>
    static Mvec<T> e0120246(){
        Mvec<T> res;
        return res.componentToOne(5, 82);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012024i2.
    template<typename T>
    static Mvec<T> e012024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 83);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0120256.
    template<typename T>
    static Mvec<T> e0120256(){
        Mvec<T> res;
        return res.componentToOne(5, 84);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012025i2.
    template<typename T>
    static Mvec<T> e012025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 85);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012026i2.
    template<typename T>
    static Mvec<T> e012026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 86);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012456.
    template<typename T>
    static Mvec<T> e012456(){
        Mvec<T> res;
        return res.componentToOne(5, 87);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01245i2.
    template<typename T>
    static Mvec<T> e01245i2(){
        Mvec<T> res;
        return res.componentToOne(5, 88);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01246i2.
    template<typename T>
    static Mvec<T> e01246i2(){
        Mvec<T> res;
        return res.componentToOne(5, 89);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01256i2.
    template<typename T>
    static Mvec<T> e01256i2(){
        Mvec<T> res;
        return res.componentToOne(5, 90);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1024.
    template<typename T>
    static Mvec<T> e013i1024(){
        Mvec<T> res;
        return res.componentToOne(5, 91);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1025.
    template<typename T>
    static Mvec<T> e013i1025(){
        Mvec<T> res;
        return res.componentToOne(5, 92);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1026.
    template<typename T>
    static Mvec<T> e013i1026(){
        Mvec<T> res;
        return res.componentToOne(5, 93);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i102i2.
    template<typename T>
    static Mvec<T> e013i102i2(){
        Mvec<T> res;
        return res.componentToOne(5, 94);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i145.
    template<typename T>
    static Mvec<T> e013i145(){
        Mvec<T> res;
        return res.componentToOne(5, 95);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i146.
    template<typename T>
    static Mvec<T> e013i146(){
        Mvec<T> res;
        return res.componentToOne(5, 96);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i14i2.
    template<typename T>
    static Mvec<T> e013i14i2(){
        Mvec<T> res;
        return res.componentToOne(5, 97);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i156.
    template<typename T>
    static Mvec<T> e013i156(){
        Mvec<T> res;
        return res.componentToOne(5, 98);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i15i2.
    template<typename T>
    static Mvec<T> e013i15i2(){
        Mvec<T> res;
        return res.componentToOne(5, 99);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i16i2.
    template<typename T>
    static Mvec<T> e013i16i2(){
        Mvec<T> res;
        return res.componentToOne(5, 100);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0130245.
    template<typename T>
    static Mvec<T> e0130245(){
        Mvec<T> res;
        return res.componentToOne(5, 101);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0130246.
    template<typename T>
    static Mvec<T> e0130246(){
        Mvec<T> res;
        return res.componentToOne(5, 102);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013024i2.
    template<typename T>
    static Mvec<T> e013024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 103);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0130256.
    template<typename T>
    static Mvec<T> e0130256(){
        Mvec<T> res;
        return res.componentToOne(5, 104);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013025i2.
    template<typename T>
    static Mvec<T> e013025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 105);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013026i2.
    template<typename T>
    static Mvec<T> e013026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 106);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013456.
    template<typename T>
    static Mvec<T> e013456(){
        Mvec<T> res;
        return res.componentToOne(5, 107);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01345i2.
    template<typename T>
    static Mvec<T> e01345i2(){
        Mvec<T> res;
        return res.componentToOne(5, 108);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01346i2.
    template<typename T>
    static Mvec<T> e01346i2(){
        Mvec<T> res;
        return res.componentToOne(5, 109);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01356i2.
    template<typename T>
    static Mvec<T> e01356i2(){
        Mvec<T> res;
        return res.componentToOne(5, 110);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i10245.
    template<typename T>
    static Mvec<T> e01i10245(){
        Mvec<T> res;
        return res.componentToOne(5, 111);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i10246.
    template<typename T>
    static Mvec<T> e01i10246(){
        Mvec<T> res;
        return res.componentToOne(5, 112);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1024i2.
    template<typename T>
    static Mvec<T> e01i1024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 113);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i10256.
    template<typename T>
    static Mvec<T> e01i10256(){
        Mvec<T> res;
        return res.componentToOne(5, 114);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1025i2.
    template<typename T>
    static Mvec<T> e01i1025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 115);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1026i2.
    template<typename T>
    static Mvec<T> e01i1026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 116);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1456.
    template<typename T>
    static Mvec<T> e01i1456(){
        Mvec<T> res;
        return res.componentToOne(5, 117);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i145i2.
    template<typename T>
    static Mvec<T> e01i145i2(){
        Mvec<T> res;
        return res.componentToOne(5, 118);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i146i2.
    template<typename T>
    static Mvec<T> e01i146i2(){
        Mvec<T> res;
        return res.componentToOne(5, 119);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i156i2.
    template<typename T>
    static Mvec<T> e01i156i2(){
        Mvec<T> res;
        return res.componentToOne(5, 120);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0102456.
    template<typename T>
    static Mvec<T> e0102456(){
        Mvec<T> res;
        return res.componentToOne(5, 121);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 010245i2.
    template<typename T>
    static Mvec<T> e010245i2(){
        Mvec<T> res;
        return res.componentToOne(5, 122);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 010246i2.
    template<typename T>
    static Mvec<T> e010246i2(){
        Mvec<T> res;
        return res.componentToOne(5, 123);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 010256i2.
    template<typename T>
    static Mvec<T> e010256i2(){
        Mvec<T> res;
        return res.componentToOne(5, 124);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01456i2.
    template<typename T>
    static Mvec<T> e01456i2(){
        Mvec<T> res;
        return res.componentToOne(5, 125);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i102.
    template<typename T>
    static Mvec<T> e123i102(){
        Mvec<T> res;
        return res.componentToOne(5, 126);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i14.
    template<typename T>
    static Mvec<T> e123i14(){
        Mvec<T> res;
        return res.componentToOne(5, 127);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i15.
    template<typename T>
    static Mvec<T> e123i15(){
        Mvec<T> res;
        return res.componentToOne(5, 128);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i16.
    template<typename T>
    static Mvec<T> e123i16(){
        Mvec<T> res;
        return res.componentToOne(5, 129);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1i2.
    template<typename T>
    static Mvec<T> e123i1i2(){
        Mvec<T> res;
        return res.componentToOne(5, 130);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123024.
    template<typename T>
    static Mvec<T> e123024(){
        Mvec<T> res;
        return res.componentToOne(5, 131);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123025.
    template<typename T>
    static Mvec<T> e123025(){
        Mvec<T> res;
        return res.componentToOne(5, 132);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123026.
    template<typename T>
    static Mvec<T> e123026(){
        Mvec<T> res;
        return res.componentToOne(5, 133);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12302i2.
    template<typename T>
    static Mvec<T> e12302i2(){
        Mvec<T> res;
        return res.componentToOne(5, 134);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12345.
    template<typename T>
    static Mvec<T> e12345(){
        Mvec<T> res;
        return res.componentToOne(5, 135);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12346.
    template<typename T>
    static Mvec<T> e12346(){
        Mvec<T> res;
        return res.componentToOne(5, 136);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1234i2.
    template<typename T>
    static Mvec<T> e1234i2(){
        Mvec<T> res;
        return res.componentToOne(5, 137);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12356.
    template<typename T>
    static Mvec<T> e12356(){
        Mvec<T> res;
        return res.componentToOne(5, 138);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1235i2.
    template<typename T>
    static Mvec<T> e1235i2(){
        Mvec<T> res;
        return res.componentToOne(5, 139);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1236i2.
    template<typename T>
    static Mvec<T> e1236i2(){
        Mvec<T> res;
        return res.componentToOne(5, 140);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1024.
    template<typename T>
    static Mvec<T> e12i1024(){
        Mvec<T> res;
        return res.componentToOne(5, 141);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1025.
    template<typename T>
    static Mvec<T> e12i1025(){
        Mvec<T> res;
        return res.componentToOne(5, 142);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1026.
    template<typename T>
    static Mvec<T> e12i1026(){
        Mvec<T> res;
        return res.componentToOne(5, 143);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i102i2.
    template<typename T>
    static Mvec<T> e12i102i2(){
        Mvec<T> res;
        return res.componentToOne(5, 144);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i145.
    template<typename T>
    static Mvec<T> e12i145(){
        Mvec<T> res;
        return res.componentToOne(5, 145);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i146.
    template<typename T>
    static Mvec<T> e12i146(){
        Mvec<T> res;
        return res.componentToOne(5, 146);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i14i2.
    template<typename T>
    static Mvec<T> e12i14i2(){
        Mvec<T> res;
        return res.componentToOne(5, 147);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i156.
    template<typename T>
    static Mvec<T> e12i156(){
        Mvec<T> res;
        return res.componentToOne(5, 148);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i15i2.
    template<typename T>
    static Mvec<T> e12i15i2(){
        Mvec<T> res;
        return res.componentToOne(5, 149);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i16i2.
    template<typename T>
    static Mvec<T> e12i16i2(){
        Mvec<T> res;
        return res.componentToOne(5, 150);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 120245.
    template<typename T>
    static Mvec<T> e120245(){
        Mvec<T> res;
        return res.componentToOne(5, 151);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 120246.
    template<typename T>
    static Mvec<T> e120246(){
        Mvec<T> res;
        return res.componentToOne(5, 152);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12024i2.
    template<typename T>
    static Mvec<T> e12024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 153);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 120256.
    template<typename T>
    static Mvec<T> e120256(){
        Mvec<T> res;
        return res.componentToOne(5, 154);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12025i2.
    template<typename T>
    static Mvec<T> e12025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 155);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12026i2.
    template<typename T>
    static Mvec<T> e12026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 156);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12456.
    template<typename T>
    static Mvec<T> e12456(){
        Mvec<T> res;
        return res.componentToOne(5, 157);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1245i2.
    template<typename T>
    static Mvec<T> e1245i2(){
        Mvec<T> res;
        return res.componentToOne(5, 158);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1246i2.
    template<typename T>
    static Mvec<T> e1246i2(){
        Mvec<T> res;
        return res.componentToOne(5, 159);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1256i2.
    template<typename T>
    static Mvec<T> e1256i2(){
        Mvec<T> res;
        return res.componentToOne(5, 160);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1024.
    template<typename T>
    static Mvec<T> e13i1024(){
        Mvec<T> res;
        return res.componentToOne(5, 161);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1025.
    template<typename T>
    static Mvec<T> e13i1025(){
        Mvec<T> res;
        return res.componentToOne(5, 162);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1026.
    template<typename T>
    static Mvec<T> e13i1026(){
        Mvec<T> res;
        return res.componentToOne(5, 163);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i102i2.
    template<typename T>
    static Mvec<T> e13i102i2(){
        Mvec<T> res;
        return res.componentToOne(5, 164);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i145.
    template<typename T>
    static Mvec<T> e13i145(){
        Mvec<T> res;
        return res.componentToOne(5, 165);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i146.
    template<typename T>
    static Mvec<T> e13i146(){
        Mvec<T> res;
        return res.componentToOne(5, 166);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i14i2.
    template<typename T>
    static Mvec<T> e13i14i2(){
        Mvec<T> res;
        return res.componentToOne(5, 167);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i156.
    template<typename T>
    static Mvec<T> e13i156(){
        Mvec<T> res;
        return res.componentToOne(5, 168);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i15i2.
    template<typename T>
    static Mvec<T> e13i15i2(){
        Mvec<T> res;
        return res.componentToOne(5, 169);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i16i2.
    template<typename T>
    static Mvec<T> e13i16i2(){
        Mvec<T> res;
        return res.componentToOne(5, 170);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 130245.
    template<typename T>
    static Mvec<T> e130245(){
        Mvec<T> res;
        return res.componentToOne(5, 171);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 130246.
    template<typename T>
    static Mvec<T> e130246(){
        Mvec<T> res;
        return res.componentToOne(5, 172);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13024i2.
    template<typename T>
    static Mvec<T> e13024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 173);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 130256.
    template<typename T>
    static Mvec<T> e130256(){
        Mvec<T> res;
        return res.componentToOne(5, 174);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13025i2.
    template<typename T>
    static Mvec<T> e13025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 175);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13026i2.
    template<typename T>
    static Mvec<T> e13026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 176);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13456.
    template<typename T>
    static Mvec<T> e13456(){
        Mvec<T> res;
        return res.componentToOne(5, 177);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1345i2.
    template<typename T>
    static Mvec<T> e1345i2(){
        Mvec<T> res;
        return res.componentToOne(5, 178);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1346i2.
    template<typename T>
    static Mvec<T> e1346i2(){
        Mvec<T> res;
        return res.componentToOne(5, 179);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1356i2.
    template<typename T>
    static Mvec<T> e1356i2(){
        Mvec<T> res;
        return res.componentToOne(5, 180);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i10245.
    template<typename T>
    static Mvec<T> e1i10245(){
        Mvec<T> res;
        return res.componentToOne(5, 181);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i10246.
    template<typename T>
    static Mvec<T> e1i10246(){
        Mvec<T> res;
        return res.componentToOne(5, 182);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1024i2.
    template<typename T>
    static Mvec<T> e1i1024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 183);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i10256.
    template<typename T>
    static Mvec<T> e1i10256(){
        Mvec<T> res;
        return res.componentToOne(5, 184);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1025i2.
    template<typename T>
    static Mvec<T> e1i1025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 185);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1026i2.
    template<typename T>
    static Mvec<T> e1i1026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 186);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1456.
    template<typename T>
    static Mvec<T> e1i1456(){
        Mvec<T> res;
        return res.componentToOne(5, 187);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i145i2.
    template<typename T>
    static Mvec<T> e1i145i2(){
        Mvec<T> res;
        return res.componentToOne(5, 188);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i146i2.
    template<typename T>
    static Mvec<T> e1i146i2(){
        Mvec<T> res;
        return res.componentToOne(5, 189);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i156i2.
    template<typename T>
    static Mvec<T> e1i156i2(){
        Mvec<T> res;
        return res.componentToOne(5, 190);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 102456.
    template<typename T>
    static Mvec<T> e102456(){
        Mvec<T> res;
        return res.componentToOne(5, 191);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 10245i2.
    template<typename T>
    static Mvec<T> e10245i2(){
        Mvec<T> res;
        return res.componentToOne(5, 192);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 10246i2.
    template<typename T>
    static Mvec<T> e10246i2(){
        Mvec<T> res;
        return res.componentToOne(5, 193);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 10256i2.
    template<typename T>
    static Mvec<T> e10256i2(){
        Mvec<T> res;
        return res.componentToOne(5, 194);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1456i2.
    template<typename T>
    static Mvec<T> e1456i2(){
        Mvec<T> res;
        return res.componentToOne(5, 195);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1024.
    template<typename T>
    static Mvec<T> e23i1024(){
        Mvec<T> res;
        return res.componentToOne(5, 196);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1025.
    template<typename T>
    static Mvec<T> e23i1025(){
        Mvec<T> res;
        return res.componentToOne(5, 197);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1026.
    template<typename T>
    static Mvec<T> e23i1026(){
        Mvec<T> res;
        return res.componentToOne(5, 198);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i102i2.
    template<typename T>
    static Mvec<T> e23i102i2(){
        Mvec<T> res;
        return res.componentToOne(5, 199);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i145.
    template<typename T>
    static Mvec<T> e23i145(){
        Mvec<T> res;
        return res.componentToOne(5, 200);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i146.
    template<typename T>
    static Mvec<T> e23i146(){
        Mvec<T> res;
        return res.componentToOne(5, 201);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i14i2.
    template<typename T>
    static Mvec<T> e23i14i2(){
        Mvec<T> res;
        return res.componentToOne(5, 202);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i156.
    template<typename T>
    static Mvec<T> e23i156(){
        Mvec<T> res;
        return res.componentToOne(5, 203);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i15i2.
    template<typename T>
    static Mvec<T> e23i15i2(){
        Mvec<T> res;
        return res.componentToOne(5, 204);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i16i2.
    template<typename T>
    static Mvec<T> e23i16i2(){
        Mvec<T> res;
        return res.componentToOne(5, 205);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 230245.
    template<typename T>
    static Mvec<T> e230245(){
        Mvec<T> res;
        return res.componentToOne(5, 206);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 230246.
    template<typename T>
    static Mvec<T> e230246(){
        Mvec<T> res;
        return res.componentToOne(5, 207);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23024i2.
    template<typename T>
    static Mvec<T> e23024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 208);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 230256.
    template<typename T>
    static Mvec<T> e230256(){
        Mvec<T> res;
        return res.componentToOne(5, 209);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23025i2.
    template<typename T>
    static Mvec<T> e23025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 210);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23026i2.
    template<typename T>
    static Mvec<T> e23026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 211);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23456.
    template<typename T>
    static Mvec<T> e23456(){
        Mvec<T> res;
        return res.componentToOne(5, 212);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2345i2.
    template<typename T>
    static Mvec<T> e2345i2(){
        Mvec<T> res;
        return res.componentToOne(5, 213);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2346i2.
    template<typename T>
    static Mvec<T> e2346i2(){
        Mvec<T> res;
        return res.componentToOne(5, 214);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2356i2.
    template<typename T>
    static Mvec<T> e2356i2(){
        Mvec<T> res;
        return res.componentToOne(5, 215);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i10245.
    template<typename T>
    static Mvec<T> e2i10245(){
        Mvec<T> res;
        return res.componentToOne(5, 216);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i10246.
    template<typename T>
    static Mvec<T> e2i10246(){
        Mvec<T> res;
        return res.componentToOne(5, 217);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1024i2.
    template<typename T>
    static Mvec<T> e2i1024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 218);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i10256.
    template<typename T>
    static Mvec<T> e2i10256(){
        Mvec<T> res;
        return res.componentToOne(5, 219);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1025i2.
    template<typename T>
    static Mvec<T> e2i1025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 220);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1026i2.
    template<typename T>
    static Mvec<T> e2i1026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 221);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1456.
    template<typename T>
    static Mvec<T> e2i1456(){
        Mvec<T> res;
        return res.componentToOne(5, 222);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i145i2.
    template<typename T>
    static Mvec<T> e2i145i2(){
        Mvec<T> res;
        return res.componentToOne(5, 223);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i146i2.
    template<typename T>
    static Mvec<T> e2i146i2(){
        Mvec<T> res;
        return res.componentToOne(5, 224);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i156i2.
    template<typename T>
    static Mvec<T> e2i156i2(){
        Mvec<T> res;
        return res.componentToOne(5, 225);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 202456.
    template<typename T>
    static Mvec<T> e202456(){
        Mvec<T> res;
        return res.componentToOne(5, 226);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 20245i2.
    template<typename T>
    static Mvec<T> e20245i2(){
        Mvec<T> res;
        return res.componentToOne(5, 227);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 20246i2.
    template<typename T>
    static Mvec<T> e20246i2(){
        Mvec<T> res;
        return res.componentToOne(5, 228);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 20256i2.
    template<typename T>
    static Mvec<T> e20256i2(){
        Mvec<T> res;
        return res.componentToOne(5, 229);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2456i2.
    template<typename T>
    static Mvec<T> e2456i2(){
        Mvec<T> res;
        return res.componentToOne(5, 230);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i10245.
    template<typename T>
    static Mvec<T> e3i10245(){
        Mvec<T> res;
        return res.componentToOne(5, 231);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i10246.
    template<typename T>
    static Mvec<T> e3i10246(){
        Mvec<T> res;
        return res.componentToOne(5, 232);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1024i2.
    template<typename T>
    static Mvec<T> e3i1024i2(){
        Mvec<T> res;
        return res.componentToOne(5, 233);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i10256.
    template<typename T>
    static Mvec<T> e3i10256(){
        Mvec<T> res;
        return res.componentToOne(5, 234);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1025i2.
    template<typename T>
    static Mvec<T> e3i1025i2(){
        Mvec<T> res;
        return res.componentToOne(5, 235);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1026i2.
    template<typename T>
    static Mvec<T> e3i1026i2(){
        Mvec<T> res;
        return res.componentToOne(5, 236);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1456.
    template<typename T>
    static Mvec<T> e3i1456(){
        Mvec<T> res;
        return res.componentToOne(5, 237);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i145i2.
    template<typename T>
    static Mvec<T> e3i145i2(){
        Mvec<T> res;
        return res.componentToOne(5, 238);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i146i2.
    template<typename T>
    static Mvec<T> e3i146i2(){
        Mvec<T> res;
        return res.componentToOne(5, 239);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i156i2.
    template<typename T>
    static Mvec<T> e3i156i2(){
        Mvec<T> res;
        return res.componentToOne(5, 240);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 302456.
    template<typename T>
    static Mvec<T> e302456(){
        Mvec<T> res;
        return res.componentToOne(5, 241);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 30245i2.
    template<typename T>
    static Mvec<T> e30245i2(){
        Mvec<T> res;
        return res.componentToOne(5, 242);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 30246i2.
    template<typename T>
    static Mvec<T> e30246i2(){
        Mvec<T> res;
        return res.componentToOne(5, 243);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 30256i2.
    template<typename T>
    static Mvec<T> e30256i2(){
        Mvec<T> res;
        return res.componentToOne(5, 244);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3456i2.
    template<typename T>
    static Mvec<T> e3456i2(){
        Mvec<T> res;
        return res.componentToOne(5, 245);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i102456.
    template<typename T>
    static Mvec<T> ei102456(){
        Mvec<T> res;
        return res.componentToOne(5, 246);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i10245i2.
    template<typename T>
    static Mvec<T> ei10245i2(){
        Mvec<T> res;
        return res.componentToOne(5, 247);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i10246i2.
    template<typename T>
    static Mvec<T> ei10246i2(){
        Mvec<T> res;
        return res.componentToOne(5, 248);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i10256i2.
    template<typename T>
    static Mvec<T> ei10256i2(){
        Mvec<T> res;
        return res.componentToOne(5, 249);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i1456i2.
    template<typename T>
    static Mvec<T> ei1456i2(){
        Mvec<T> res;
        return res.componentToOne(5, 250);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 02456i2.
    template<typename T>
    static Mvec<T> e02456i2(){
        Mvec<T> res;
        return res.componentToOne(5, 251);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i102.
    template<typename T>
    static Mvec<T> e01123i102(){
        Mvec<T> res;
        return res.componentToOne(6, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i14.
    template<typename T>
    static Mvec<T> e01123i14(){
        Mvec<T> res;
        return res.componentToOne(6, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i15.
    template<typename T>
    static Mvec<T> e01123i15(){
        Mvec<T> res;
        return res.componentToOne(6, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i16.
    template<typename T>
    static Mvec<T> e01123i16(){
        Mvec<T> res;
        return res.componentToOne(6, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1i2.
    template<typename T>
    static Mvec<T> e01123i1i2(){
        Mvec<T> res;
        return res.componentToOne(6, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123024.
    template<typename T>
    static Mvec<T> e01123024(){
        Mvec<T> res;
        return res.componentToOne(6, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123025.
    template<typename T>
    static Mvec<T> e01123025(){
        Mvec<T> res;
        return res.componentToOne(6, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123026.
    template<typename T>
    static Mvec<T> e01123026(){
        Mvec<T> res;
        return res.componentToOne(6, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112302i2.
    template<typename T>
    static Mvec<T> e0112302i2(){
        Mvec<T> res;
        return res.componentToOne(6, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112345.
    template<typename T>
    static Mvec<T> e0112345(){
        Mvec<T> res;
        return res.componentToOne(6, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112346.
    template<typename T>
    static Mvec<T> e0112346(){
        Mvec<T> res;
        return res.componentToOne(6, 10);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011234i2.
    template<typename T>
    static Mvec<T> e011234i2(){
        Mvec<T> res;
        return res.componentToOne(6, 11);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112356.
    template<typename T>
    static Mvec<T> e0112356(){
        Mvec<T> res;
        return res.componentToOne(6, 12);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011235i2.
    template<typename T>
    static Mvec<T> e011235i2(){
        Mvec<T> res;
        return res.componentToOne(6, 13);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011236i2.
    template<typename T>
    static Mvec<T> e011236i2(){
        Mvec<T> res;
        return res.componentToOne(6, 14);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1024.
    template<typename T>
    static Mvec<T> e0112i1024(){
        Mvec<T> res;
        return res.componentToOne(6, 15);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1025.
    template<typename T>
    static Mvec<T> e0112i1025(){
        Mvec<T> res;
        return res.componentToOne(6, 16);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1026.
    template<typename T>
    static Mvec<T> e0112i1026(){
        Mvec<T> res;
        return res.componentToOne(6, 17);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i102i2.
    template<typename T>
    static Mvec<T> e0112i102i2(){
        Mvec<T> res;
        return res.componentToOne(6, 18);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i145.
    template<typename T>
    static Mvec<T> e0112i145(){
        Mvec<T> res;
        return res.componentToOne(6, 19);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i146.
    template<typename T>
    static Mvec<T> e0112i146(){
        Mvec<T> res;
        return res.componentToOne(6, 20);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i14i2.
    template<typename T>
    static Mvec<T> e0112i14i2(){
        Mvec<T> res;
        return res.componentToOne(6, 21);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i156.
    template<typename T>
    static Mvec<T> e0112i156(){
        Mvec<T> res;
        return res.componentToOne(6, 22);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i15i2.
    template<typename T>
    static Mvec<T> e0112i15i2(){
        Mvec<T> res;
        return res.componentToOne(6, 23);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i16i2.
    template<typename T>
    static Mvec<T> e0112i16i2(){
        Mvec<T> res;
        return res.componentToOne(6, 24);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01120245.
    template<typename T>
    static Mvec<T> e01120245(){
        Mvec<T> res;
        return res.componentToOne(6, 25);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01120246.
    template<typename T>
    static Mvec<T> e01120246(){
        Mvec<T> res;
        return res.componentToOne(6, 26);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112024i2.
    template<typename T>
    static Mvec<T> e0112024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 27);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01120256.
    template<typename T>
    static Mvec<T> e01120256(){
        Mvec<T> res;
        return res.componentToOne(6, 28);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112025i2.
    template<typename T>
    static Mvec<T> e0112025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 29);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112026i2.
    template<typename T>
    static Mvec<T> e0112026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 30);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112456.
    template<typename T>
    static Mvec<T> e0112456(){
        Mvec<T> res;
        return res.componentToOne(6, 31);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011245i2.
    template<typename T>
    static Mvec<T> e011245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 32);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011246i2.
    template<typename T>
    static Mvec<T> e011246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 33);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011256i2.
    template<typename T>
    static Mvec<T> e011256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 34);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1024.
    template<typename T>
    static Mvec<T> e0113i1024(){
        Mvec<T> res;
        return res.componentToOne(6, 35);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1025.
    template<typename T>
    static Mvec<T> e0113i1025(){
        Mvec<T> res;
        return res.componentToOne(6, 36);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1026.
    template<typename T>
    static Mvec<T> e0113i1026(){
        Mvec<T> res;
        return res.componentToOne(6, 37);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i102i2.
    template<typename T>
    static Mvec<T> e0113i102i2(){
        Mvec<T> res;
        return res.componentToOne(6, 38);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i145.
    template<typename T>
    static Mvec<T> e0113i145(){
        Mvec<T> res;
        return res.componentToOne(6, 39);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i146.
    template<typename T>
    static Mvec<T> e0113i146(){
        Mvec<T> res;
        return res.componentToOne(6, 40);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i14i2.
    template<typename T>
    static Mvec<T> e0113i14i2(){
        Mvec<T> res;
        return res.componentToOne(6, 41);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i156.
    template<typename T>
    static Mvec<T> e0113i156(){
        Mvec<T> res;
        return res.componentToOne(6, 42);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i15i2.
    template<typename T>
    static Mvec<T> e0113i15i2(){
        Mvec<T> res;
        return res.componentToOne(6, 43);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i16i2.
    template<typename T>
    static Mvec<T> e0113i16i2(){
        Mvec<T> res;
        return res.componentToOne(6, 44);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01130245.
    template<typename T>
    static Mvec<T> e01130245(){
        Mvec<T> res;
        return res.componentToOne(6, 45);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01130246.
    template<typename T>
    static Mvec<T> e01130246(){
        Mvec<T> res;
        return res.componentToOne(6, 46);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113024i2.
    template<typename T>
    static Mvec<T> e0113024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 47);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01130256.
    template<typename T>
    static Mvec<T> e01130256(){
        Mvec<T> res;
        return res.componentToOne(6, 48);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113025i2.
    template<typename T>
    static Mvec<T> e0113025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 49);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113026i2.
    template<typename T>
    static Mvec<T> e0113026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 50);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113456.
    template<typename T>
    static Mvec<T> e0113456(){
        Mvec<T> res;
        return res.componentToOne(6, 51);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011345i2.
    template<typename T>
    static Mvec<T> e011345i2(){
        Mvec<T> res;
        return res.componentToOne(6, 52);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011346i2.
    template<typename T>
    static Mvec<T> e011346i2(){
        Mvec<T> res;
        return res.componentToOne(6, 53);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011356i2.
    template<typename T>
    static Mvec<T> e011356i2(){
        Mvec<T> res;
        return res.componentToOne(6, 54);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i10245.
    template<typename T>
    static Mvec<T> e011i10245(){
        Mvec<T> res;
        return res.componentToOne(6, 55);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i10246.
    template<typename T>
    static Mvec<T> e011i10246(){
        Mvec<T> res;
        return res.componentToOne(6, 56);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1024i2.
    template<typename T>
    static Mvec<T> e011i1024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 57);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i10256.
    template<typename T>
    static Mvec<T> e011i10256(){
        Mvec<T> res;
        return res.componentToOne(6, 58);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1025i2.
    template<typename T>
    static Mvec<T> e011i1025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 59);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1026i2.
    template<typename T>
    static Mvec<T> e011i1026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 60);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1456.
    template<typename T>
    static Mvec<T> e011i1456(){
        Mvec<T> res;
        return res.componentToOne(6, 61);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i145i2.
    template<typename T>
    static Mvec<T> e011i145i2(){
        Mvec<T> res;
        return res.componentToOne(6, 62);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i146i2.
    template<typename T>
    static Mvec<T> e011i146i2(){
        Mvec<T> res;
        return res.componentToOne(6, 63);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i156i2.
    template<typename T>
    static Mvec<T> e011i156i2(){
        Mvec<T> res;
        return res.componentToOne(6, 64);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01102456.
    template<typename T>
    static Mvec<T> e01102456(){
        Mvec<T> res;
        return res.componentToOne(6, 65);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0110245i2.
    template<typename T>
    static Mvec<T> e0110245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 66);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0110246i2.
    template<typename T>
    static Mvec<T> e0110246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 67);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0110256i2.
    template<typename T>
    static Mvec<T> e0110256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 68);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011456i2.
    template<typename T>
    static Mvec<T> e011456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 69);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1024.
    template<typename T>
    static Mvec<T> e0123i1024(){
        Mvec<T> res;
        return res.componentToOne(6, 70);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1025.
    template<typename T>
    static Mvec<T> e0123i1025(){
        Mvec<T> res;
        return res.componentToOne(6, 71);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1026.
    template<typename T>
    static Mvec<T> e0123i1026(){
        Mvec<T> res;
        return res.componentToOne(6, 72);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i102i2.
    template<typename T>
    static Mvec<T> e0123i102i2(){
        Mvec<T> res;
        return res.componentToOne(6, 73);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i145.
    template<typename T>
    static Mvec<T> e0123i145(){
        Mvec<T> res;
        return res.componentToOne(6, 74);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i146.
    template<typename T>
    static Mvec<T> e0123i146(){
        Mvec<T> res;
        return res.componentToOne(6, 75);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i14i2.
    template<typename T>
    static Mvec<T> e0123i14i2(){
        Mvec<T> res;
        return res.componentToOne(6, 76);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i156.
    template<typename T>
    static Mvec<T> e0123i156(){
        Mvec<T> res;
        return res.componentToOne(6, 77);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i15i2.
    template<typename T>
    static Mvec<T> e0123i15i2(){
        Mvec<T> res;
        return res.componentToOne(6, 78);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i16i2.
    template<typename T>
    static Mvec<T> e0123i16i2(){
        Mvec<T> res;
        return res.componentToOne(6, 79);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01230245.
    template<typename T>
    static Mvec<T> e01230245(){
        Mvec<T> res;
        return res.componentToOne(6, 80);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01230246.
    template<typename T>
    static Mvec<T> e01230246(){
        Mvec<T> res;
        return res.componentToOne(6, 81);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123024i2.
    template<typename T>
    static Mvec<T> e0123024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 82);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01230256.
    template<typename T>
    static Mvec<T> e01230256(){
        Mvec<T> res;
        return res.componentToOne(6, 83);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123025i2.
    template<typename T>
    static Mvec<T> e0123025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 84);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123026i2.
    template<typename T>
    static Mvec<T> e0123026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 85);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123456.
    template<typename T>
    static Mvec<T> e0123456(){
        Mvec<T> res;
        return res.componentToOne(6, 86);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012345i2.
    template<typename T>
    static Mvec<T> e012345i2(){
        Mvec<T> res;
        return res.componentToOne(6, 87);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012346i2.
    template<typename T>
    static Mvec<T> e012346i2(){
        Mvec<T> res;
        return res.componentToOne(6, 88);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012356i2.
    template<typename T>
    static Mvec<T> e012356i2(){
        Mvec<T> res;
        return res.componentToOne(6, 89);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i10245.
    template<typename T>
    static Mvec<T> e012i10245(){
        Mvec<T> res;
        return res.componentToOne(6, 90);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i10246.
    template<typename T>
    static Mvec<T> e012i10246(){
        Mvec<T> res;
        return res.componentToOne(6, 91);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1024i2.
    template<typename T>
    static Mvec<T> e012i1024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 92);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i10256.
    template<typename T>
    static Mvec<T> e012i10256(){
        Mvec<T> res;
        return res.componentToOne(6, 93);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1025i2.
    template<typename T>
    static Mvec<T> e012i1025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 94);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1026i2.
    template<typename T>
    static Mvec<T> e012i1026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 95);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1456.
    template<typename T>
    static Mvec<T> e012i1456(){
        Mvec<T> res;
        return res.componentToOne(6, 96);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i145i2.
    template<typename T>
    static Mvec<T> e012i145i2(){
        Mvec<T> res;
        return res.componentToOne(6, 97);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i146i2.
    template<typename T>
    static Mvec<T> e012i146i2(){
        Mvec<T> res;
        return res.componentToOne(6, 98);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i156i2.
    template<typename T>
    static Mvec<T> e012i156i2(){
        Mvec<T> res;
        return res.componentToOne(6, 99);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01202456.
    template<typename T>
    static Mvec<T> e01202456(){
        Mvec<T> res;
        return res.componentToOne(6, 100);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0120245i2.
    template<typename T>
    static Mvec<T> e0120245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 101);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0120246i2.
    template<typename T>
    static Mvec<T> e0120246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 102);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0120256i2.
    template<typename T>
    static Mvec<T> e0120256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 103);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012456i2.
    template<typename T>
    static Mvec<T> e012456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 104);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i10245.
    template<typename T>
    static Mvec<T> e013i10245(){
        Mvec<T> res;
        return res.componentToOne(6, 105);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i10246.
    template<typename T>
    static Mvec<T> e013i10246(){
        Mvec<T> res;
        return res.componentToOne(6, 106);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1024i2.
    template<typename T>
    static Mvec<T> e013i1024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 107);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i10256.
    template<typename T>
    static Mvec<T> e013i10256(){
        Mvec<T> res;
        return res.componentToOne(6, 108);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1025i2.
    template<typename T>
    static Mvec<T> e013i1025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 109);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1026i2.
    template<typename T>
    static Mvec<T> e013i1026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 110);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1456.
    template<typename T>
    static Mvec<T> e013i1456(){
        Mvec<T> res;
        return res.componentToOne(6, 111);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i145i2.
    template<typename T>
    static Mvec<T> e013i145i2(){
        Mvec<T> res;
        return res.componentToOne(6, 112);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i146i2.
    template<typename T>
    static Mvec<T> e013i146i2(){
        Mvec<T> res;
        return res.componentToOne(6, 113);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i156i2.
    template<typename T>
    static Mvec<T> e013i156i2(){
        Mvec<T> res;
        return res.componentToOne(6, 114);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01302456.
    template<typename T>
    static Mvec<T> e01302456(){
        Mvec<T> res;
        return res.componentToOne(6, 115);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0130245i2.
    template<typename T>
    static Mvec<T> e0130245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 116);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0130246i2.
    template<typename T>
    static Mvec<T> e0130246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 117);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0130256i2.
    template<typename T>
    static Mvec<T> e0130256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 118);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013456i2.
    template<typename T>
    static Mvec<T> e013456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 119);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i102456.
    template<typename T>
    static Mvec<T> e01i102456(){
        Mvec<T> res;
        return res.componentToOne(6, 120);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i10245i2.
    template<typename T>
    static Mvec<T> e01i10245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 121);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i10246i2.
    template<typename T>
    static Mvec<T> e01i10246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 122);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i10256i2.
    template<typename T>
    static Mvec<T> e01i10256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 123);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i1456i2.
    template<typename T>
    static Mvec<T> e01i1456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 124);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0102456i2.
    template<typename T>
    static Mvec<T> e0102456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 125);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1024.
    template<typename T>
    static Mvec<T> e123i1024(){
        Mvec<T> res;
        return res.componentToOne(6, 126);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1025.
    template<typename T>
    static Mvec<T> e123i1025(){
        Mvec<T> res;
        return res.componentToOne(6, 127);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1026.
    template<typename T>
    static Mvec<T> e123i1026(){
        Mvec<T> res;
        return res.componentToOne(6, 128);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i102i2.
    template<typename T>
    static Mvec<T> e123i102i2(){
        Mvec<T> res;
        return res.componentToOne(6, 129);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i145.
    template<typename T>
    static Mvec<T> e123i145(){
        Mvec<T> res;
        return res.componentToOne(6, 130);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i146.
    template<typename T>
    static Mvec<T> e123i146(){
        Mvec<T> res;
        return res.componentToOne(6, 131);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i14i2.
    template<typename T>
    static Mvec<T> e123i14i2(){
        Mvec<T> res;
        return res.componentToOne(6, 132);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i156.
    template<typename T>
    static Mvec<T> e123i156(){
        Mvec<T> res;
        return res.componentToOne(6, 133);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i15i2.
    template<typename T>
    static Mvec<T> e123i15i2(){
        Mvec<T> res;
        return res.componentToOne(6, 134);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i16i2.
    template<typename T>
    static Mvec<T> e123i16i2(){
        Mvec<T> res;
        return res.componentToOne(6, 135);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1230245.
    template<typename T>
    static Mvec<T> e1230245(){
        Mvec<T> res;
        return res.componentToOne(6, 136);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1230246.
    template<typename T>
    static Mvec<T> e1230246(){
        Mvec<T> res;
        return res.componentToOne(6, 137);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123024i2.
    template<typename T>
    static Mvec<T> e123024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 138);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1230256.
    template<typename T>
    static Mvec<T> e1230256(){
        Mvec<T> res;
        return res.componentToOne(6, 139);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123025i2.
    template<typename T>
    static Mvec<T> e123025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 140);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123026i2.
    template<typename T>
    static Mvec<T> e123026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 141);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123456.
    template<typename T>
    static Mvec<T> e123456(){
        Mvec<T> res;
        return res.componentToOne(6, 142);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12345i2.
    template<typename T>
    static Mvec<T> e12345i2(){
        Mvec<T> res;
        return res.componentToOne(6, 143);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12346i2.
    template<typename T>
    static Mvec<T> e12346i2(){
        Mvec<T> res;
        return res.componentToOne(6, 144);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12356i2.
    template<typename T>
    static Mvec<T> e12356i2(){
        Mvec<T> res;
        return res.componentToOne(6, 145);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i10245.
    template<typename T>
    static Mvec<T> e12i10245(){
        Mvec<T> res;
        return res.componentToOne(6, 146);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i10246.
    template<typename T>
    static Mvec<T> e12i10246(){
        Mvec<T> res;
        return res.componentToOne(6, 147);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1024i2.
    template<typename T>
    static Mvec<T> e12i1024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 148);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i10256.
    template<typename T>
    static Mvec<T> e12i10256(){
        Mvec<T> res;
        return res.componentToOne(6, 149);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1025i2.
    template<typename T>
    static Mvec<T> e12i1025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 150);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1026i2.
    template<typename T>
    static Mvec<T> e12i1026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 151);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1456.
    template<typename T>
    static Mvec<T> e12i1456(){
        Mvec<T> res;
        return res.componentToOne(6, 152);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i145i2.
    template<typename T>
    static Mvec<T> e12i145i2(){
        Mvec<T> res;
        return res.componentToOne(6, 153);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i146i2.
    template<typename T>
    static Mvec<T> e12i146i2(){
        Mvec<T> res;
        return res.componentToOne(6, 154);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i156i2.
    template<typename T>
    static Mvec<T> e12i156i2(){
        Mvec<T> res;
        return res.componentToOne(6, 155);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1202456.
    template<typename T>
    static Mvec<T> e1202456(){
        Mvec<T> res;
        return res.componentToOne(6, 156);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 120245i2.
    template<typename T>
    static Mvec<T> e120245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 157);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 120246i2.
    template<typename T>
    static Mvec<T> e120246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 158);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 120256i2.
    template<typename T>
    static Mvec<T> e120256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 159);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12456i2.
    template<typename T>
    static Mvec<T> e12456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 160);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i10245.
    template<typename T>
    static Mvec<T> e13i10245(){
        Mvec<T> res;
        return res.componentToOne(6, 161);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i10246.
    template<typename T>
    static Mvec<T> e13i10246(){
        Mvec<T> res;
        return res.componentToOne(6, 162);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1024i2.
    template<typename T>
    static Mvec<T> e13i1024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 163);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i10256.
    template<typename T>
    static Mvec<T> e13i10256(){
        Mvec<T> res;
        return res.componentToOne(6, 164);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1025i2.
    template<typename T>
    static Mvec<T> e13i1025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 165);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1026i2.
    template<typename T>
    static Mvec<T> e13i1026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 166);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1456.
    template<typename T>
    static Mvec<T> e13i1456(){
        Mvec<T> res;
        return res.componentToOne(6, 167);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i145i2.
    template<typename T>
    static Mvec<T> e13i145i2(){
        Mvec<T> res;
        return res.componentToOne(6, 168);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i146i2.
    template<typename T>
    static Mvec<T> e13i146i2(){
        Mvec<T> res;
        return res.componentToOne(6, 169);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i156i2.
    template<typename T>
    static Mvec<T> e13i156i2(){
        Mvec<T> res;
        return res.componentToOne(6, 170);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1302456.
    template<typename T>
    static Mvec<T> e1302456(){
        Mvec<T> res;
        return res.componentToOne(6, 171);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 130245i2.
    template<typename T>
    static Mvec<T> e130245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 172);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 130246i2.
    template<typename T>
    static Mvec<T> e130246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 173);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 130256i2.
    template<typename T>
    static Mvec<T> e130256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 174);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13456i2.
    template<typename T>
    static Mvec<T> e13456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 175);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i102456.
    template<typename T>
    static Mvec<T> e1i102456(){
        Mvec<T> res;
        return res.componentToOne(6, 176);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i10245i2.
    template<typename T>
    static Mvec<T> e1i10245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 177);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i10246i2.
    template<typename T>
    static Mvec<T> e1i10246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 178);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i10256i2.
    template<typename T>
    static Mvec<T> e1i10256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 179);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i1456i2.
    template<typename T>
    static Mvec<T> e1i1456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 180);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 102456i2.
    template<typename T>
    static Mvec<T> e102456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 181);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i10245.
    template<typename T>
    static Mvec<T> e23i10245(){
        Mvec<T> res;
        return res.componentToOne(6, 182);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i10246.
    template<typename T>
    static Mvec<T> e23i10246(){
        Mvec<T> res;
        return res.componentToOne(6, 183);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1024i2.
    template<typename T>
    static Mvec<T> e23i1024i2(){
        Mvec<T> res;
        return res.componentToOne(6, 184);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i10256.
    template<typename T>
    static Mvec<T> e23i10256(){
        Mvec<T> res;
        return res.componentToOne(6, 185);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1025i2.
    template<typename T>
    static Mvec<T> e23i1025i2(){
        Mvec<T> res;
        return res.componentToOne(6, 186);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1026i2.
    template<typename T>
    static Mvec<T> e23i1026i2(){
        Mvec<T> res;
        return res.componentToOne(6, 187);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1456.
    template<typename T>
    static Mvec<T> e23i1456(){
        Mvec<T> res;
        return res.componentToOne(6, 188);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i145i2.
    template<typename T>
    static Mvec<T> e23i145i2(){
        Mvec<T> res;
        return res.componentToOne(6, 189);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i146i2.
    template<typename T>
    static Mvec<T> e23i146i2(){
        Mvec<T> res;
        return res.componentToOne(6, 190);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i156i2.
    template<typename T>
    static Mvec<T> e23i156i2(){
        Mvec<T> res;
        return res.componentToOne(6, 191);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2302456.
    template<typename T>
    static Mvec<T> e2302456(){
        Mvec<T> res;
        return res.componentToOne(6, 192);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 230245i2.
    template<typename T>
    static Mvec<T> e230245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 193);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 230246i2.
    template<typename T>
    static Mvec<T> e230246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 194);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 230256i2.
    template<typename T>
    static Mvec<T> e230256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 195);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23456i2.
    template<typename T>
    static Mvec<T> e23456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 196);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i102456.
    template<typename T>
    static Mvec<T> e2i102456(){
        Mvec<T> res;
        return res.componentToOne(6, 197);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i10245i2.
    template<typename T>
    static Mvec<T> e2i10245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 198);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i10246i2.
    template<typename T>
    static Mvec<T> e2i10246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 199);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i10256i2.
    template<typename T>
    static Mvec<T> e2i10256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 200);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i1456i2.
    template<typename T>
    static Mvec<T> e2i1456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 201);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 202456i2.
    template<typename T>
    static Mvec<T> e202456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 202);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i102456.
    template<typename T>
    static Mvec<T> e3i102456(){
        Mvec<T> res;
        return res.componentToOne(6, 203);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i10245i2.
    template<typename T>
    static Mvec<T> e3i10245i2(){
        Mvec<T> res;
        return res.componentToOne(6, 204);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i10246i2.
    template<typename T>
    static Mvec<T> e3i10246i2(){
        Mvec<T> res;
        return res.componentToOne(6, 205);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i10256i2.
    template<typename T>
    static Mvec<T> e3i10256i2(){
        Mvec<T> res;
        return res.componentToOne(6, 206);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i1456i2.
    template<typename T>
    static Mvec<T> e3i1456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 207);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 302456i2.
    template<typename T>
    static Mvec<T> e302456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 208);
    }

    /// \brief return a multivector that contains only the unit basis k-vector i102456i2.
    template<typename T>
    static Mvec<T> ei102456i2(){
        Mvec<T> res;
        return res.componentToOne(6, 209);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1024.
    template<typename T>
    static Mvec<T> e01123i1024(){
        Mvec<T> res;
        return res.componentToOne(7, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1025.
    template<typename T>
    static Mvec<T> e01123i1025(){
        Mvec<T> res;
        return res.componentToOne(7, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1026.
    template<typename T>
    static Mvec<T> e01123i1026(){
        Mvec<T> res;
        return res.componentToOne(7, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i102i2.
    template<typename T>
    static Mvec<T> e01123i102i2(){
        Mvec<T> res;
        return res.componentToOne(7, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i145.
    template<typename T>
    static Mvec<T> e01123i145(){
        Mvec<T> res;
        return res.componentToOne(7, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i146.
    template<typename T>
    static Mvec<T> e01123i146(){
        Mvec<T> res;
        return res.componentToOne(7, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i14i2.
    template<typename T>
    static Mvec<T> e01123i14i2(){
        Mvec<T> res;
        return res.componentToOne(7, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i156.
    template<typename T>
    static Mvec<T> e01123i156(){
        Mvec<T> res;
        return res.componentToOne(7, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i15i2.
    template<typename T>
    static Mvec<T> e01123i15i2(){
        Mvec<T> res;
        return res.componentToOne(7, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i16i2.
    template<typename T>
    static Mvec<T> e01123i16i2(){
        Mvec<T> res;
        return res.componentToOne(7, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011230245.
    template<typename T>
    static Mvec<T> e011230245(){
        Mvec<T> res;
        return res.componentToOne(7, 10);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011230246.
    template<typename T>
    static Mvec<T> e011230246(){
        Mvec<T> res;
        return res.componentToOne(7, 11);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123024i2.
    template<typename T>
    static Mvec<T> e01123024i2(){
        Mvec<T> res;
        return res.componentToOne(7, 12);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011230256.
    template<typename T>
    static Mvec<T> e011230256(){
        Mvec<T> res;
        return res.componentToOne(7, 13);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123025i2.
    template<typename T>
    static Mvec<T> e01123025i2(){
        Mvec<T> res;
        return res.componentToOne(7, 14);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123026i2.
    template<typename T>
    static Mvec<T> e01123026i2(){
        Mvec<T> res;
        return res.componentToOne(7, 15);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123456.
    template<typename T>
    static Mvec<T> e01123456(){
        Mvec<T> res;
        return res.componentToOne(7, 16);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112345i2.
    template<typename T>
    static Mvec<T> e0112345i2(){
        Mvec<T> res;
        return res.componentToOne(7, 17);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112346i2.
    template<typename T>
    static Mvec<T> e0112346i2(){
        Mvec<T> res;
        return res.componentToOne(7, 18);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112356i2.
    template<typename T>
    static Mvec<T> e0112356i2(){
        Mvec<T> res;
        return res.componentToOne(7, 19);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i10245.
    template<typename T>
    static Mvec<T> e0112i10245(){
        Mvec<T> res;
        return res.componentToOne(7, 20);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i10246.
    template<typename T>
    static Mvec<T> e0112i10246(){
        Mvec<T> res;
        return res.componentToOne(7, 21);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1024i2.
    template<typename T>
    static Mvec<T> e0112i1024i2(){
        Mvec<T> res;
        return res.componentToOne(7, 22);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i10256.
    template<typename T>
    static Mvec<T> e0112i10256(){
        Mvec<T> res;
        return res.componentToOne(7, 23);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1025i2.
    template<typename T>
    static Mvec<T> e0112i1025i2(){
        Mvec<T> res;
        return res.componentToOne(7, 24);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1026i2.
    template<typename T>
    static Mvec<T> e0112i1026i2(){
        Mvec<T> res;
        return res.componentToOne(7, 25);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1456.
    template<typename T>
    static Mvec<T> e0112i1456(){
        Mvec<T> res;
        return res.componentToOne(7, 26);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i145i2.
    template<typename T>
    static Mvec<T> e0112i145i2(){
        Mvec<T> res;
        return res.componentToOne(7, 27);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i146i2.
    template<typename T>
    static Mvec<T> e0112i146i2(){
        Mvec<T> res;
        return res.componentToOne(7, 28);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i156i2.
    template<typename T>
    static Mvec<T> e0112i156i2(){
        Mvec<T> res;
        return res.componentToOne(7, 29);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011202456.
    template<typename T>
    static Mvec<T> e011202456(){
        Mvec<T> res;
        return res.componentToOne(7, 30);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01120245i2.
    template<typename T>
    static Mvec<T> e01120245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 31);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01120246i2.
    template<typename T>
    static Mvec<T> e01120246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 32);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01120256i2.
    template<typename T>
    static Mvec<T> e01120256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 33);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112456i2.
    template<typename T>
    static Mvec<T> e0112456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 34);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i10245.
    template<typename T>
    static Mvec<T> e0113i10245(){
        Mvec<T> res;
        return res.componentToOne(7, 35);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i10246.
    template<typename T>
    static Mvec<T> e0113i10246(){
        Mvec<T> res;
        return res.componentToOne(7, 36);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1024i2.
    template<typename T>
    static Mvec<T> e0113i1024i2(){
        Mvec<T> res;
        return res.componentToOne(7, 37);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i10256.
    template<typename T>
    static Mvec<T> e0113i10256(){
        Mvec<T> res;
        return res.componentToOne(7, 38);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1025i2.
    template<typename T>
    static Mvec<T> e0113i1025i2(){
        Mvec<T> res;
        return res.componentToOne(7, 39);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1026i2.
    template<typename T>
    static Mvec<T> e0113i1026i2(){
        Mvec<T> res;
        return res.componentToOne(7, 40);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1456.
    template<typename T>
    static Mvec<T> e0113i1456(){
        Mvec<T> res;
        return res.componentToOne(7, 41);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i145i2.
    template<typename T>
    static Mvec<T> e0113i145i2(){
        Mvec<T> res;
        return res.componentToOne(7, 42);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i146i2.
    template<typename T>
    static Mvec<T> e0113i146i2(){
        Mvec<T> res;
        return res.componentToOne(7, 43);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i156i2.
    template<typename T>
    static Mvec<T> e0113i156i2(){
        Mvec<T> res;
        return res.componentToOne(7, 44);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011302456.
    template<typename T>
    static Mvec<T> e011302456(){
        Mvec<T> res;
        return res.componentToOne(7, 45);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01130245i2.
    template<typename T>
    static Mvec<T> e01130245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 46);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01130246i2.
    template<typename T>
    static Mvec<T> e01130246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 47);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01130256i2.
    template<typename T>
    static Mvec<T> e01130256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 48);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113456i2.
    template<typename T>
    static Mvec<T> e0113456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 49);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i102456.
    template<typename T>
    static Mvec<T> e011i102456(){
        Mvec<T> res;
        return res.componentToOne(7, 50);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i10245i2.
    template<typename T>
    static Mvec<T> e011i10245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 51);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i10246i2.
    template<typename T>
    static Mvec<T> e011i10246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 52);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i10256i2.
    template<typename T>
    static Mvec<T> e011i10256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 53);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i1456i2.
    template<typename T>
    static Mvec<T> e011i1456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 54);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01102456i2.
    template<typename T>
    static Mvec<T> e01102456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 55);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i10245.
    template<typename T>
    static Mvec<T> e0123i10245(){
        Mvec<T> res;
        return res.componentToOne(7, 56);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i10246.
    template<typename T>
    static Mvec<T> e0123i10246(){
        Mvec<T> res;
        return res.componentToOne(7, 57);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1024i2.
    template<typename T>
    static Mvec<T> e0123i1024i2(){
        Mvec<T> res;
        return res.componentToOne(7, 58);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i10256.
    template<typename T>
    static Mvec<T> e0123i10256(){
        Mvec<T> res;
        return res.componentToOne(7, 59);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1025i2.
    template<typename T>
    static Mvec<T> e0123i1025i2(){
        Mvec<T> res;
        return res.componentToOne(7, 60);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1026i2.
    template<typename T>
    static Mvec<T> e0123i1026i2(){
        Mvec<T> res;
        return res.componentToOne(7, 61);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1456.
    template<typename T>
    static Mvec<T> e0123i1456(){
        Mvec<T> res;
        return res.componentToOne(7, 62);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i145i2.
    template<typename T>
    static Mvec<T> e0123i145i2(){
        Mvec<T> res;
        return res.componentToOne(7, 63);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i146i2.
    template<typename T>
    static Mvec<T> e0123i146i2(){
        Mvec<T> res;
        return res.componentToOne(7, 64);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i156i2.
    template<typename T>
    static Mvec<T> e0123i156i2(){
        Mvec<T> res;
        return res.componentToOne(7, 65);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012302456.
    template<typename T>
    static Mvec<T> e012302456(){
        Mvec<T> res;
        return res.componentToOne(7, 66);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01230245i2.
    template<typename T>
    static Mvec<T> e01230245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 67);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01230246i2.
    template<typename T>
    static Mvec<T> e01230246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 68);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01230256i2.
    template<typename T>
    static Mvec<T> e01230256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 69);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123456i2.
    template<typename T>
    static Mvec<T> e0123456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 70);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i102456.
    template<typename T>
    static Mvec<T> e012i102456(){
        Mvec<T> res;
        return res.componentToOne(7, 71);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i10245i2.
    template<typename T>
    static Mvec<T> e012i10245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 72);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i10246i2.
    template<typename T>
    static Mvec<T> e012i10246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 73);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i10256i2.
    template<typename T>
    static Mvec<T> e012i10256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 74);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i1456i2.
    template<typename T>
    static Mvec<T> e012i1456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 75);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01202456i2.
    template<typename T>
    static Mvec<T> e01202456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 76);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i102456.
    template<typename T>
    static Mvec<T> e013i102456(){
        Mvec<T> res;
        return res.componentToOne(7, 77);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i10245i2.
    template<typename T>
    static Mvec<T> e013i10245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 78);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i10246i2.
    template<typename T>
    static Mvec<T> e013i10246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 79);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i10256i2.
    template<typename T>
    static Mvec<T> e013i10256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 80);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i1456i2.
    template<typename T>
    static Mvec<T> e013i1456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 81);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01302456i2.
    template<typename T>
    static Mvec<T> e01302456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 82);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01i102456i2.
    template<typename T>
    static Mvec<T> e01i102456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 83);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i10245.
    template<typename T>
    static Mvec<T> e123i10245(){
        Mvec<T> res;
        return res.componentToOne(7, 84);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i10246.
    template<typename T>
    static Mvec<T> e123i10246(){
        Mvec<T> res;
        return res.componentToOne(7, 85);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1024i2.
    template<typename T>
    static Mvec<T> e123i1024i2(){
        Mvec<T> res;
        return res.componentToOne(7, 86);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i10256.
    template<typename T>
    static Mvec<T> e123i10256(){
        Mvec<T> res;
        return res.componentToOne(7, 87);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1025i2.
    template<typename T>
    static Mvec<T> e123i1025i2(){
        Mvec<T> res;
        return res.componentToOne(7, 88);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1026i2.
    template<typename T>
    static Mvec<T> e123i1026i2(){
        Mvec<T> res;
        return res.componentToOne(7, 89);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1456.
    template<typename T>
    static Mvec<T> e123i1456(){
        Mvec<T> res;
        return res.componentToOne(7, 90);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i145i2.
    template<typename T>
    static Mvec<T> e123i145i2(){
        Mvec<T> res;
        return res.componentToOne(7, 91);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i146i2.
    template<typename T>
    static Mvec<T> e123i146i2(){
        Mvec<T> res;
        return res.componentToOne(7, 92);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i156i2.
    template<typename T>
    static Mvec<T> e123i156i2(){
        Mvec<T> res;
        return res.componentToOne(7, 93);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12302456.
    template<typename T>
    static Mvec<T> e12302456(){
        Mvec<T> res;
        return res.componentToOne(7, 94);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1230245i2.
    template<typename T>
    static Mvec<T> e1230245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 95);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1230246i2.
    template<typename T>
    static Mvec<T> e1230246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 96);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1230256i2.
    template<typename T>
    static Mvec<T> e1230256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 97);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123456i2.
    template<typename T>
    static Mvec<T> e123456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 98);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i102456.
    template<typename T>
    static Mvec<T> e12i102456(){
        Mvec<T> res;
        return res.componentToOne(7, 99);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i10245i2.
    template<typename T>
    static Mvec<T> e12i10245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 100);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i10246i2.
    template<typename T>
    static Mvec<T> e12i10246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 101);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i10256i2.
    template<typename T>
    static Mvec<T> e12i10256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 102);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i1456i2.
    template<typename T>
    static Mvec<T> e12i1456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 103);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1202456i2.
    template<typename T>
    static Mvec<T> e1202456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 104);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i102456.
    template<typename T>
    static Mvec<T> e13i102456(){
        Mvec<T> res;
        return res.componentToOne(7, 105);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i10245i2.
    template<typename T>
    static Mvec<T> e13i10245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 106);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i10246i2.
    template<typename T>
    static Mvec<T> e13i10246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 107);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i10256i2.
    template<typename T>
    static Mvec<T> e13i10256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 108);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i1456i2.
    template<typename T>
    static Mvec<T> e13i1456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 109);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1302456i2.
    template<typename T>
    static Mvec<T> e1302456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 110);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 1i102456i2.
    template<typename T>
    static Mvec<T> e1i102456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 111);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i102456.
    template<typename T>
    static Mvec<T> e23i102456(){
        Mvec<T> res;
        return res.componentToOne(7, 112);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i10245i2.
    template<typename T>
    static Mvec<T> e23i10245i2(){
        Mvec<T> res;
        return res.componentToOne(7, 113);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i10246i2.
    template<typename T>
    static Mvec<T> e23i10246i2(){
        Mvec<T> res;
        return res.componentToOne(7, 114);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i10256i2.
    template<typename T>
    static Mvec<T> e23i10256i2(){
        Mvec<T> res;
        return res.componentToOne(7, 115);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i1456i2.
    template<typename T>
    static Mvec<T> e23i1456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 116);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2302456i2.
    template<typename T>
    static Mvec<T> e2302456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 117);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 2i102456i2.
    template<typename T>
    static Mvec<T> e2i102456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 118);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 3i102456i2.
    template<typename T>
    static Mvec<T> e3i102456i2(){
        Mvec<T> res;
        return res.componentToOne(7, 119);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i10245.
    template<typename T>
    static Mvec<T> e01123i10245(){
        Mvec<T> res;
        return res.componentToOne(8, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i10246.
    template<typename T>
    static Mvec<T> e01123i10246(){
        Mvec<T> res;
        return res.componentToOne(8, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1024i2.
    template<typename T>
    static Mvec<T> e01123i1024i2(){
        Mvec<T> res;
        return res.componentToOne(8, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i10256.
    template<typename T>
    static Mvec<T> e01123i10256(){
        Mvec<T> res;
        return res.componentToOne(8, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1025i2.
    template<typename T>
    static Mvec<T> e01123i1025i2(){
        Mvec<T> res;
        return res.componentToOne(8, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1026i2.
    template<typename T>
    static Mvec<T> e01123i1026i2(){
        Mvec<T> res;
        return res.componentToOne(8, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1456.
    template<typename T>
    static Mvec<T> e01123i1456(){
        Mvec<T> res;
        return res.componentToOne(8, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i145i2.
    template<typename T>
    static Mvec<T> e01123i145i2(){
        Mvec<T> res;
        return res.componentToOne(8, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i146i2.
    template<typename T>
    static Mvec<T> e01123i146i2(){
        Mvec<T> res;
        return res.componentToOne(8, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i156i2.
    template<typename T>
    static Mvec<T> e01123i156i2(){
        Mvec<T> res;
        return res.componentToOne(8, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112302456.
    template<typename T>
    static Mvec<T> e0112302456(){
        Mvec<T> res;
        return res.componentToOne(8, 10);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011230245i2.
    template<typename T>
    static Mvec<T> e011230245i2(){
        Mvec<T> res;
        return res.componentToOne(8, 11);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011230246i2.
    template<typename T>
    static Mvec<T> e011230246i2(){
        Mvec<T> res;
        return res.componentToOne(8, 12);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011230256i2.
    template<typename T>
    static Mvec<T> e011230256i2(){
        Mvec<T> res;
        return res.componentToOne(8, 13);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123456i2.
    template<typename T>
    static Mvec<T> e01123456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 14);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i102456.
    template<typename T>
    static Mvec<T> e0112i102456(){
        Mvec<T> res;
        return res.componentToOne(8, 15);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i10245i2.
    template<typename T>
    static Mvec<T> e0112i10245i2(){
        Mvec<T> res;
        return res.componentToOne(8, 16);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i10246i2.
    template<typename T>
    static Mvec<T> e0112i10246i2(){
        Mvec<T> res;
        return res.componentToOne(8, 17);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i10256i2.
    template<typename T>
    static Mvec<T> e0112i10256i2(){
        Mvec<T> res;
        return res.componentToOne(8, 18);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i1456i2.
    template<typename T>
    static Mvec<T> e0112i1456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 19);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011202456i2.
    template<typename T>
    static Mvec<T> e011202456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 20);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i102456.
    template<typename T>
    static Mvec<T> e0113i102456(){
        Mvec<T> res;
        return res.componentToOne(8, 21);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i10245i2.
    template<typename T>
    static Mvec<T> e0113i10245i2(){
        Mvec<T> res;
        return res.componentToOne(8, 22);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i10246i2.
    template<typename T>
    static Mvec<T> e0113i10246i2(){
        Mvec<T> res;
        return res.componentToOne(8, 23);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i10256i2.
    template<typename T>
    static Mvec<T> e0113i10256i2(){
        Mvec<T> res;
        return res.componentToOne(8, 24);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i1456i2.
    template<typename T>
    static Mvec<T> e0113i1456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 25);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011302456i2.
    template<typename T>
    static Mvec<T> e011302456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 26);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 011i102456i2.
    template<typename T>
    static Mvec<T> e011i102456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 27);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i102456.
    template<typename T>
    static Mvec<T> e0123i102456(){
        Mvec<T> res;
        return res.componentToOne(8, 28);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i10245i2.
    template<typename T>
    static Mvec<T> e0123i10245i2(){
        Mvec<T> res;
        return res.componentToOne(8, 29);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i10246i2.
    template<typename T>
    static Mvec<T> e0123i10246i2(){
        Mvec<T> res;
        return res.componentToOne(8, 30);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i10256i2.
    template<typename T>
    static Mvec<T> e0123i10256i2(){
        Mvec<T> res;
        return res.componentToOne(8, 31);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i1456i2.
    template<typename T>
    static Mvec<T> e0123i1456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 32);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012302456i2.
    template<typename T>
    static Mvec<T> e012302456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 33);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 012i102456i2.
    template<typename T>
    static Mvec<T> e012i102456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 34);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 013i102456i2.
    template<typename T>
    static Mvec<T> e013i102456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 35);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i102456.
    template<typename T>
    static Mvec<T> e123i102456(){
        Mvec<T> res;
        return res.componentToOne(8, 36);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i10245i2.
    template<typename T>
    static Mvec<T> e123i10245i2(){
        Mvec<T> res;
        return res.componentToOne(8, 37);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i10246i2.
    template<typename T>
    static Mvec<T> e123i10246i2(){
        Mvec<T> res;
        return res.componentToOne(8, 38);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i10256i2.
    template<typename T>
    static Mvec<T> e123i10256i2(){
        Mvec<T> res;
        return res.componentToOne(8, 39);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i1456i2.
    template<typename T>
    static Mvec<T> e123i1456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 40);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12302456i2.
    template<typename T>
    static Mvec<T> e12302456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 41);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 12i102456i2.
    template<typename T>
    static Mvec<T> e12i102456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 42);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 13i102456i2.
    template<typename T>
    static Mvec<T> e13i102456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 43);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 23i102456i2.
    template<typename T>
    static Mvec<T> e23i102456i2(){
        Mvec<T> res;
        return res.componentToOne(8, 44);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i102456.
    template<typename T>
    static Mvec<T> e01123i102456(){
        Mvec<T> res;
        return res.componentToOne(9, 0);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i10245i2.
    template<typename T>
    static Mvec<T> e01123i10245i2(){
        Mvec<T> res;
        return res.componentToOne(9, 1);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i10246i2.
    template<typename T>
    static Mvec<T> e01123i10246i2(){
        Mvec<T> res;
        return res.componentToOne(9, 2);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i10256i2.
    template<typename T>
    static Mvec<T> e01123i10256i2(){
        Mvec<T> res;
        return res.componentToOne(9, 3);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i1456i2.
    template<typename T>
    static Mvec<T> e01123i1456i2(){
        Mvec<T> res;
        return res.componentToOne(9, 4);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112302456i2.
    template<typename T>
    static Mvec<T> e0112302456i2(){
        Mvec<T> res;
        return res.componentToOne(9, 5);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0112i102456i2.
    template<typename T>
    static Mvec<T> e0112i102456i2(){
        Mvec<T> res;
        return res.componentToOne(9, 6);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0113i102456i2.
    template<typename T>
    static Mvec<T> e0113i102456i2(){
        Mvec<T> res;
        return res.componentToOne(9, 7);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 0123i102456i2.
    template<typename T>
    static Mvec<T> e0123i102456i2(){
        Mvec<T> res;
        return res.componentToOne(9, 8);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 123i102456i2.
    template<typename T>
    static Mvec<T> e123i102456i2(){
        Mvec<T> res;
        return res.componentToOne(9, 9);
    }

    /// \brief return a multivector that contains only the unit basis k-vector 01123i102456i2.
    template<typename T>
    static Mvec<T> e01123i102456i2(){
        Mvec<T> res;
        return res.componentToOne(10, 0);
    }



    //    template<typename U>
    //    void recursiveTraversalMultivector(std::ostream &stream, const Mvec<U> &mvec, unsigned int currentGrade, int currentIndex, std::vector<int> listBasisBlades, unsigned int lastIndex, unsigned int gradeMV, bool& moreThanOne);


    void temporaryFunction1();

}     /// End of Namespace

#endif // C3GA2_MULTI_VECTOR_HPP__
