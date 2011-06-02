#ifndef __itkIntermodesThresholdImageFilter_txx
#define __itkIntermodesThresholdImageFilter_txx
#include "itkIntermodesThresholdImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkIntermodesThresholdImageCalculator.h"
#include "itkProgressAccumulator.h"

namespace itk {

template<class TInputImage, class TOutputImage>
IntermodesThresholdImageFilter<TInputImage, TOutputImage>
::IntermodesThresholdImageFilter()
{
  m_OutsideValue   = NumericTraits<OutputPixelType>::Zero;
  m_InsideValue    = NumericTraits<OutputPixelType>::max();
  m_Threshold      = NumericTraits<InputPixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_MaxSmoothingIterations = 10000;
}

template<class TInputImage, class TOutputImage>
void
IntermodesThresholdImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  typename ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  // Compute the Intermodes Threshold for the input image
  typename IntermodesThresholdImageCalculator<TInputImage>::Pointer calculator =
    IntermodesThresholdImageCalculator<TInputImage>::New();
  calculator->SetImage (this->GetInput());
  calculator->SetNumberOfHistogramBins (m_NumberOfHistogramBins);
  calculator->SetMaxSmoothingIterations(m_MaxSmoothingIterations);
  calculator->Compute();
  m_Threshold = calculator->GetThreshold();

  typename BinaryThresholdImageFilter<TInputImage,TOutputImage>::Pointer threshold =
    BinaryThresholdImageFilter<TInputImage,TOutputImage>::New();

  progress->RegisterInternalFilter(threshold,.5f);
  threshold->GraftOutput (this->GetOutput());
  threshold->SetInput (this->GetInput());
  threshold->SetLowerThreshold(NumericTraits<InputPixelType>::NonpositiveMin());
  threshold->SetUpperThreshold(calculator->GetThreshold());
  threshold->SetInsideValue (m_InsideValue);
  threshold->SetOutsideValue (m_OutsideValue);
  threshold->Update();

  this->GraftOutput(threshold->GetOutput());
}

template<class TInputImage, class TOutputImage>
void
IntermodesThresholdImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  TInputImage * input = const_cast<TInputImage *>(this->GetInput());
  if( input )
    {
    input->SetRequestedRegionToLargestPossibleRegion();
    }
}

template<class TInputImage, class TOutputImage>
void 
IntermodesThresholdImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "OutsideValue: "
     << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_OutsideValue) << std::endl;
  os << indent << "InsideValue: "
     << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_InsideValue) << std::endl;
  os << indent << "NumberOfHistogramBins: "
     << m_NumberOfHistogramBins << std::endl;
  os << indent << "Threshold (computed): "
     << static_cast<typename NumericTraits<InputPixelType>::PrintType>(m_Threshold) << std::endl;

}


}// end namespace itk
#endif
