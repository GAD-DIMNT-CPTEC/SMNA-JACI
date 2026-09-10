
=========================================================================================

Instrucoes para instalar, compilar e executar o Modelo Global no SX6, na Una ou ClusterPC

=========================================================================================

-------------------------
1 Estrutura de diretorios
-------------------------


                                         ....
                                           |
                                           |
                  ------------------------------------------------------
                  |                                                    |
                model                                                 pos
                  |                                                    |
       -----------------------                              ------------------------
       |       |      |       |                             |        |      |      |
    datain  dataout  exec  GLOBALMPI                      datain  dataout  exec  POSTGRIB
               |                                                     |  
       ----------------                                       ---------------- 
       |              |                                       |              |    
    T042L28  ....  T341L64                                 T042L28 .....  T341L64
    
    
  1.1 O diretorio "model/datain" armazena os dados necessarios para uma simulacao.
  1.2 O diretorio "model/dataout/<resolucao>" onde <resolucao> eh na forma T042L28 
      armazena os arquivos de saida.
  1.3 O diretorio "model/exec" armazena o executavel do modelo, o script de execucao
      e o namelist MODELIN. O MODELIN deve ter os "paths" completos dos diretorios
      model/datain e model/dataout.
  1.4 O diretorio "model/GLOBALMPI" armazena os fontes do modelo global.
  1.5 Os diretorios para o pos-processamento seguem as mesmas regras.

  1.6 Crie os diretorios "model" com os sub-diretorios "datain", "dataout" e "exec.
      O diretorio GLOBALMPI serah criado automaticamente durante a obtencao do fonte.
  1.7 Crie o sub-diretorio "dataout/<resolucao>" para a <resolucao> desejada.
  1.8 Repita o procedimento para o pos.


=========================================================================================

---------------------------------
2 Obtencao dos fontes do modelo
---------------------------------
    
  2.1 Dentro do diretorio model, obtenha os fontes do cvs:
      "  cvs checkout GLOBALMPI "
      o subdiretorio GLOBALMPI eh gerado automaticamente.

=========================================================================================

-----------------------------------------
3 Compilacao dos fontes do Modelo Global
-----------------------------------------

  3.0 Dentro do sub-diretorio GLOBALMPI, edite o arquivo Makefile. Retire o primeiro caracter
      (#) das linhas relativas a maquina desejada. Por exemplo, para o SX-6, transforme as linhas
      -----------------------------------------------------------------------------------
#
#  SX6
#
#FTRACE= 
#OPENMP=
#F90=sxmpif90 $(FTRACE) $(OPENMP) 
#NOASSUME=-Wf"-pvctl noassume vwork=stack"
#INLINE_CU_GRELL=-pi exp=es5
#INLINE_CU_KUO=-pi exp=es
#INLINE_CU_RAS=-pi exp=es3
#INLINE_SOUZA=-pi exp=es2
#LOADFLAG= -Wl"-Z 4G" 
      -----------------------------------------------------------------------------------
      nas linhas
      ----------------------------------------------------------------------------------
#
#  SX6
#
FTRACE= 
OPENMP=
F90=sxmpif90 $(FTRACE) $(OPENMP) 
NOASSUME=-Wf"-pvctl noassume vwork=stack"
INLINE_CU_GRELL=-pi exp=es5
INLINE_CU_KUO=-pi exp=es
INLINE_CU_RAS=-pi exp=es3
INLINE_SOUZA=-pi exp=es2
LOADFLAG= -Wl"-Z 4G" 
      ----------------------------------------------------------------------------------

  3.1 Dentro do diretorio GLOBALMPI, execute o comando:
        " make "
        Resulta arquivo executavel ParModel_MPI, no mesmo diretorio.
                                                   

=========================================================================================

--------------------------------
4 Como executar o modelo
--------------------------------

  4.1 Coloque os arquivos necessarios no diretorio model/exec. Situado no diretorio model/exec,
      execute (observe o ponto ao fim do comando):
      " ln -s ../GLOBALMPI/ParModel_MPI . "
      " cp ../GLOBALMPI/run_multi_<maquina> . "
      " cp ../GLOBALMPI/MODELIN . "
      onde <maquina> eh um de SX6, TX7 ou UNA

  4.2 Coloque os arquivos independentes de resolucao no diretorio model/datain:
      aunits, cnf2tb, cnftbl, desirtable, looktb, sibalb, sibveg, vunits

  4.3 Coloque os arquivos dependentes de resolucao no diretorio model/datain.
      Por exemplo, para a resolucao T021L09:
      GANLNMC.........unf.T021L09
      NMI.T021L09
      orgvar.T021
      sibmsk.T021
      snowfd............unf.T021
      soilms.T021
      sstaoi.T021 ou sstwkl.....T021
      t3zrl.T021
      onde .... representa o intervalo de datas de integracao desejado.

  4.4 Edite a copia do arquivo MODELIN residente no diretorio model/exec. Ajuste
      as resolucoes horizontal e vertical, o dt de integracao, as datas de
      inicio e  termino da simulacao e os diretorios de dados de entrada e
      de saida.

  4.5 Do diretorio model/exec, execute o script run_multi_<maquina>, utilizando
      o seguinte procedimento:

      (a) Caso SX6:
      Execute o comando
      " run_multi_SX6 "

      (b) Caso UNA:
      Execute o comando
      " run_multi_UNA cpu_mpi cpu_node nome hold "
      onde:
       cpu_mpi: inteiro: o numero de processos MPI desejado;
       cpu_node: inteiro: o numero de processadores por noh do cluster a utilizar
       nome: character, o nome do job (para Sun Grid Engine apenas)
       hold: argumento opcional. Se nao fornecido, o script retorna apos submeter
             o job na fila. Se fornecido, o script retorna apenas apos o termino
             do job submetido.
      Por exemplo, " run_multi_UNA 10 2 T021L09 1 " utiliza 10 processadores, alocando
      2 processadores por noh (totalizando 5 nos), submetendo o job a fila adequada e
      retornando o controle apenas apos o termino do job.

      (c) Caso TX7:
      Execute o comando
      " run_multi_TX7 cpu_mpi "
      onde:
       cpu_mpi: inteiro: o numero de processos MPI desejado;

      









