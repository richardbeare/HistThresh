#ifndef __itkIsoDataThresholdImageFilter_txx
#define __itkIsoDataThresholdImageFilter_txx
#include "itkIsoDataThresholdImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkIsoDataThresholdImageCalculator.h"
#include "itkProgressAccumulator.h"

namespace itk {

template<class TInputImage, class TOutputImage>
IsoDataThresholdImageFilter<TInputImage, TOutputImage>
::IsoDataThresholdImageFilter()
{
  m_OutsideValue   = NumericTraits<OutputPixelType>::Zero;
  m_InsideValue    = NumericTraits<OutputPixelType>::max();
  m_Threshold      = NumericTraits<InputPixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_LowThresh = 0.01;
  m_HighThresh= 0.99;
}

template<class TInputImage, class TOutputImage>
void
IsoDataThresholdImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  typename ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  // Compute the IsoData Threshold for the input image
  typename IsoDataThresholdImageCalculator<TInputImage>::Pointer triang =
    IsoDataThresholdImageCalculator<TInputImage>::New();
  triang->SetImage (this->GetInput());
  triang->SetNumberOfHistogramBins (m_NumberOfHistogramBins);
  triang->SetLowThresh(m_LowThresh);
  triang->SetHighThresh(m_HighThresh);
  triang->Compute();
  m_Threshold = triang->GetThreshold();

  typename BinaryThresholdImageFilter<TInputImage,TOutputImage>::Pointer threshold =
    BinaryThresholdImageFilter<TInputImage,TOutputImage>::New();

  progress->RegisterInternalFilter(threshold,.5f);
  threshold->GraftOutput (this->GetOutput());
  threshold->SetInput (this->GetInput());
  threshold->SetLowerThreshold(NumericTraits<InputPixelType>::NonpositiveMin());
  threshold->SetUpperThreshold(triang->GetThreshold());
  threshold->SetInsideValue (m_InsideValue);
  threshold->SetOutsideValue (m_OutsideValue);
  threshold->Update();

  this->GraftOutput(threshold->GetOutput());
}

template<class TInputImage, class TOutputImage>
void
IsoDataThresholdImageFilter<TInputImage, TOutputImage>
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
IsoDataThresholdImageFilter<TInputImage,TOutputImage>
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
