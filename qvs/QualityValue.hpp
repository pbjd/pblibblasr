#ifndef _BLASR_QUALITY_VALUE_HPP_
#define _BLASR_QUALITY_VALUE_HPP_


typedef unsigned char QualityValue;
typedef float QualityProbability;

#define MIN_QUALITY_VALUE 0
#define MAX_QUALITY_VALUE 255

enum QVScale {POverOneMinusP, // popularized by Illumina
              PHRED};

QualityValue ProbabilityToQualityValue(
    QualityProbability pErr, QVScale qvScale=POverOneMinusP); 

QualityValue PacBioQVToPhred(QualityValue pbQV); 

QualityValue ToPhred(QualityValue qv, QVScale qvScale=POverOneMinusP); 

QualityProbability QualityValueToProbability(QualityValue qv, 
    QVScale qvScale=POverOneMinusP); 
	
//static QVScale DetermineQVScaleFromChangeListID(ChangeListID &cl); 

#endif // _BLASR_QUALITY_VALUE_HPP_
