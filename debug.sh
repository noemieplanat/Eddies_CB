#!/bin/bash
#SBATCH --account=ctb-cdufour
#SBATCH --time 10:0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

for d in {14,20,27,32,29}
do
    for m in {1..12..1}
    do
        case $m in
            [0-9])
                m="0$m"
                ;;
            *)
        esac
        mv Detected_CREG_1d_OWthr_SENSI_0.3_m_122000_${m}_depthi${d}REF12_3_eddies.pickle Detected_CREG_1d_OWthr_SENSI_0.3_m_12_2000_${m}_depthi${d}REF12_3_eddies.pickle
        mv Detected_CREG_1d_OWthr_SENSI_0.3_m_242000_${m}_depthi${d}REF12_3_eddies.pickle Detected_CREG_1d_OWthr_SENSI_0.3_m_24_2000_${m}_depthi${d}REF12_3_eddies.pickle
        mv Detected_CREG_1d_OWthr_SENSI_0.3_m_742000_${m}_depthi${d}REF12_3_eddies.pickle Detected_CREG_1d_OWthr_SENSI_0.3_m_74_2000_${m}_depthi${d}REF12_3_eddies.pickle
        mv Detected_CREG_1d_OWthr_SENSI_0.3_m_1002000_${m}_depthi${d}REF12_3_eddies.pickle Detected_CREG_1d_OWthr_SENSI_0.3_m_100_2000_${m}_depthi${d}REF12_3_eddies.pickle
    done
done
    