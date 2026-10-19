/*!
        @file    channel_impl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy:  $

        @date    $LastChangedDate: #$

        @version $LastChangedRevision:  $
*/

#include "channel.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC.h"

template<typename ALLOCATOR>
void Channel_impl<ALLOCATOR>::tidyup()
{
  if( m_ptr  ) {
    BridgeACC::exit_data_delete(m_ptr, m_size);
  }
  m_ptr = nullptr;
  m_size = 0;
}


template<typename ALLOCATOR>
int Channel_impl<ALLOCATOR>::send_init(int count, int idir, int ipm)
{
#ifdef USE_ACCEL_CUDA
  m_buf.resize(1);
  bool host_has_data=false;
#else
  m_buf.resize(count);
  bool host_has_data=true;
#endif
  m_ptr = &m_buf[0];
  m_size = count;
  BridgeACC::enter_data_create(m_ptr, m_size, host_has_data);

  return m_comm.send_init(count, idir, ipm, (void *)m_ptr);
}


template<typename ALLOCATOR>
int Channel_impl<ALLOCATOR>::recv_init(int count, int idir, int ipm)
{
#ifdef USE_ACCEL_CUDA
  m_buf.resize(1);
  bool host_has_data=false;
#else
  m_buf.resize(count);
  bool host_has_data=true;
#endif
  m_ptr = &m_buf[0];
  m_size = count;
  BridgeACC::enter_data_create(m_ptr, m_size, host_has_data);
  return m_comm.recv_init(count, idir, ipm, (void *)m_ptr);
}

template class Channel_impl<std::allocator<double> >;
template class Channel_impl<std::allocator<float> >;
template class Channel_impl<std::allocator<char> >;
