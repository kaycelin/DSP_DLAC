# DSP_DLAC
1. Generate DLAC 
  input: DL signal      
  - DLconfig.bw_Channel = '20MHz'
  - DLconfig.fs = 122.88e6; % flag_fs_extend      
  ![image](https://user-images.githubusercontent.com/87049112/140466486-b8f90395-8f9a-44bc-a0bd-9eba8e7cf4c0.png)
  
  input: AC signal      
  - fsAC = [23.04e6, 122.88e6]
  - df = 100
  - Nbr = 4
  - dACSubcarrier = 2X0e3;
  - bwChannelAC = {'20MHz',fsAC,dACSubcarrier}
  ![image](https://user-images.githubusercontent.com/87049112/140466783-6f5e015d-8c4f-458f-bd32-5502499cb6bc.png)

  input: channel filter
  
  input: DLAC output power
  
  input: fNCO for [DL1 AC DL2]
  - fNCO = [-20e6 0 20e6]
  - flag_SUMCarrier = 'SUM'

  Export: DLAC waveform        
  ![image](https://user-images.githubusercontent.com/87049112/140467132-7db1085c-0894-4f4e-a2a2-7e3af5455a6e.png)

2. Branch phase shifter
  input: Branch phase shifter 
  - phsShiftNbrDeg = [0 10 -7 2]

3. Antenna combination for all branch combination to one branch

4. AWGN
  input:AWGN 
  - SNRdB_FullBand = 50     
  ![image](https://user-images.githubusercontent.com/87049112/140468456-a6cd9211-21c3-4f0e-9e11-77ffd8fa2e98.png)
  
5. Downconversion by LO (ignored)

6. NCO shift to Zero and multicarrier divide

7. DLAC channel filter for demodulation     
![image](https://user-images.githubusercontent.com/87049112/140468520-c1ba6c72-3edb-4017-9534-664685e982bc.png)

8. AC Demodulation      
![image](https://user-images.githubusercontent.com/87049112/140468560-94701978-d295-4bda-aba8-a65d142f07e6.png)

9. AC Calculate EQ coefficents and Phase Alignment
  input: AC EQ number of taps
  - ACEQ_numFirTaps = 2
  - ACEQ_NoOfAlignBr = 1
  ![image](https://user-images.githubusercontent.com/87049112/140468659-8105bf96-7b70-458a-a96a-f674f84bff3b.png)
  
10. DL(OFDM) demodulation     
![image](https://user-images.githubusercontent.com/87049112/140468785-c3480c49-4822-40ff-ac6a-e0d9c328c833.png)

11. Analysis
  - SNRin vs SNRout     
  ![image](https://user-images.githubusercontent.com/87049112/140468894-ef59f3ce-e9ed-4b4b-b57a-9b940d172ad8.png)

  - SNRin vs phaseDriftdeg      
  ![image](https://user-images.githubusercontent.com/87049112/140468950-adb7d7c0-3bbf-4075-9adf-f8fdac906b09.png)

  - SNRin vs phaseShifterForEQ      
  ![image](https://user-images.githubusercontent.com/87049112/140469046-0355eba4-34be-4bb4-a59c-d8538cb9b0cb.png)

  - SNRin vs EVM(DL)      
  ![image](https://user-images.githubusercontent.com/87049112/140469389-b01946fc-1a79-4651-a86c-984e0fd0c9b9.png) 
