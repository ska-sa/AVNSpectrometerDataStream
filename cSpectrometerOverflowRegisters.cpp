#include "cSpectrometerOverflowRegisters.h"

cSpectrometerOverflowRegisters::cSpectrometerOverflowRegisters(uint32_t u32RoachOverflowRegisters)
{
    m_bADC0OverRange = ( (u32RoachOverflowRegisters & 0x00000001) == 0x00000001 );
    m_bADC1OverRange = ( (u32RoachOverflowRegisters & 0x00000002) == 0x00000002 );
    m_bCoarseFFT0Overflow = ( (u32RoachOverflowRegisters & 0x00000004) == 0x00000004 );
    m_bCoarseFFT1Overflow = ( (u32RoachOverflowRegisters & 0x00000008) == 0x00000008 );
    m_bPacketiserOverflow = ( (u32RoachOverflowRegisters & 0x00000010) == 0x00000010 );
    m_bQueueFull = ( (u32RoachOverflowRegisters & 0x00000020) == 0x00000020 );
    m_b10GbETxOverflow = ( (u32RoachOverflowRegisters & 0x00000040) == 0x00000040 );
}

