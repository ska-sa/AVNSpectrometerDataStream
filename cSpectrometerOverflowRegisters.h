#ifndef SPECTROMETER_OVERFLOW_REGISTERS_H
#define SPECTROMETER_OVERFLOW_REGISTERS_H

#include <inttypes.h>

struct cSpectrometerOverflowRegisters
{
    cSpectrometerOverflowRegisters(uint32_t u32RoachOverflowRegisters);

    bool m_bADC0OverRange;
    bool m_bADC1OverRange;
    bool m_bCoarseFFT0Overflow;
    bool m_bCoarseFFT1Overflow;
    bool m_bPacketiserOverflow;
    bool m_bQueueFull;
    bool m_b10GbETxOverflow;
};

#endif // SPECTROMETER_OVERFLOW_REGISTERS_H
