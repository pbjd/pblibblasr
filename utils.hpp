#ifndef _BLASR_UTILS_HPP_
#define _BLASR_UTILS_HPP_

template<typename t_file>
void CrucialOpen(std::string &fileName, t_file &file, 
    std::ios_base::openmode mode=(std::ios_base::openmode)0 ); 

template<typename T_Int>
T_Int CeilOfFraction(T_Int num, T_Int denom);

template<typename T>
T*ProtectedNew(unsigned long size); 

#endif // _BLASR_UTILS_HPP_
