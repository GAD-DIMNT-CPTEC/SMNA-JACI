# Versão SMNA da branch SVN ajustada no GitHub

Repositório de origem no SVN: SMNA_v3.0.0.t12717:  
https://svn.cptec.inpe.br/smna/branch/SMNA_v3.0.0.t12717

Versão instalada na Egeon que está sendo ajustada para a JACI.

Seguir os passos:

1. Observe os pré requisitos antes de iniciar

   GitHub:  
   Pré-requisito: Git LFS  
   Caso não tenha instalado fazer --->  
   Linux (Ubuntu/Debian): `sudo apt install git-lfs`  
   macOS (Homebrew): `brew install git-lfs`  
   Windows: Baixar o instalador direto do site oficial.  
   ATIVAR: `git lfs install`

2. Clone a versão depois de atender os requisitos acima:
 
   ```
   cd /p/projetos/monan_das/${USER};
   git clone https://github.com/viezelc/SMNA_v3.t12717.git SMNA_v3.0.0.t12717;
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717;
   ```

3. Apontar para a branch de desenvolvimento da versão para a JACI:
   
   ```
   git checkout RepoSMNA-JACI;
   ```

4. Depois do repositório clonado fazer o lfs pull (passo importante):
   ```
   git lfs pull
   ```
   Obs.: se o comando acima não funcionar verifique se o git lfs está ativo no repositório SMNA_v3.0.0.t12717.
   Para ativar faça
   ```
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717;
   git lfs install;
   git lfs pull
   ```
 
5. Configuração do SMNA:
   ```
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG;
   ./config_smg.ksh configure
   ```
**Obs1.** Caso tenha ou queira mais de uma versão edite arquivo `jaci_paths.conf` e ajustar a variável "nome_smg" para um outro nome desejado. A atual versão está como "SMNA_v3.0.0.t12717/SMG"

6. Compilação do GSI e BAM:
   ```
   remova de seu login os modules load setados para que a lista seja a mais proxima da original da maquina
   vim .bashrc.jaci
 
   e remova a chamada do conda e comente os modules load
 
   comente
   o module purge || true

   veja como em https://github.com/viezelc/SMNA_v3.t12717/issues/18
 
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG;
   ./config_smg.ksh compile 
   ```

7. Verifique se todos os executavies estão presentes: Ver a lista abaixo na sequencia, Pre Bam e Pos, GSI e inctime.
   ```
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG;
   ls cptec/bam/pre/build/ParPre_MPI;
   ls cptec/bam/model/build/ParModel_MPI;
   ls cptec/bam/pos/source/POSTIN-GRIB
   ls cptec/bin/gsi.x
   ls cptec/bin/inctime
   ```

8. Testcase para o caso de completar a compilação:
   ```
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG;
   ./config_smg.ksh testcase
   ```
   Escolher opção [2].

9. Execução do pré na rodada anterior para preparação do ciclo de assimilação:
   ```
   cd  /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/cptec/bam/run;
   ./runPre -t 299 -l 64 -I 2025050900 -n 0 -O -T -G -Gt Netcdf -s
   ```

   OBS. Verificar se os arquivos necessarios serão corretamente encontrados para essa data: 2025050900. Caso dê erro por falta de arquivos uma copia esta no diretorio abaixo.
   Copia para seu /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/datainout/bam/pre/datain/ 
   ```
   ls /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/datainout/bam/pre/datain/
   ```
   
10. Rodar o Modelo para essa data anterior para preparar os FirstGuess do inicio do ciclo de assimilação:
   ```
   ./runModel -t 299 -l 64 -I 2025050900 -F 2025050909 -ts 3 -py SMT -px CPT -das -r
   ```

11. Testar o ciclo de assimilação no SMNA  
   ```
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/run;
   chmod 755 ./run_cycle.sh;
   ./run_cycle.sh -t 299 -l 64 -gt 299 -p CPT -I 2025050906 -F 2025050912
   ```
12. Os Arquivos de saí­da:
    
    ARQUIVOS setout modelo ->  
    modelo: /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/datainout/bam/model/DAS
    
    ARQUIVOS datarun GSI ->  
    GSI: /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/datarun/gsi
    
    RESULTADOS ->  
    Saídas do GSI: /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/datainout/gsi/dataout  
    Saí­das do modelo: /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/datainout/bam/model/dataout/TQ0299L064/DAS

13. Pós-processamento das previsões (rodar o Pós):  
    Obs.: verificar se o arquivo `POSTIN-GRIB.template` encontra-se no diretório `SMG/cptec/bam/run`. Este arquivo foi adicionado recentemente e pode ser encontrado em https://projetos.cptec.inpe.br/projects/smna/repository/revisions/162/entry/branch/SMNA_v3.0.x/SMG/cptec/bam/run/POSTIN-GRIB.template

   Exemplo para previsÃ£o de 5 dias: 
   ```
   cd /p/projetos/monan_das/${USER}/SMNA_v3.0.0.t12717/SMG/cptec/bam/run
   ```
Previsões:
   ``` 
   ./runModel -t 299 -l 64 -I 2025050906 -F 2025051406 -ts 6 -py CPT -px CPT
   ```
Pós:
   ```  
   ./runPos -t 299 -l 64 -I 2025050906 -F 2025051406
   ```

