# Resultados experimentais

Esta pasta reúne os dados e logs usados nas tabelas e gráficos do relatório.

## Dados consolidados

- [Baseline serial](dados/serial_raw.tsv)
- [Strong scaling](dados/strong_raw.tsv)
- [Weak scaling](dados/weak_raw.tsv)
- [Tabela calculada do strong scaling](dados/tabela_strong_calculada.tsv)
- [Tabela calculada do weak scaling](dados/tabela_weak_calculada.tsv)

## Logs por experimento

- [`logs/serial`](logs/serial): duas repetições da baseline serial e suas medições externas.
- [`logs/strong`](logs/strong): P=1, 2, 3, 4, 8, 12 e 13, com duas repetições e medições externas.
- [`logs/weak`](logs/weak): P=1, 2, 4, 8 e 13, com uma repetição e medição externa.
- [`logs/otimizacoes`](logs/otimizacoes): versão não bloqueante, paleta com LUT e regiões com MPI-IO.
- [`logs/hashes`](logs/hashes): comparações SHA-256 das versões otimizadas contra a serial.
- [`logs/execucao`](logs/execucao): logs gerais dos scripts e execuções em `nohup`.

Os arquivos com sufixo `_externo.txt` foram produzidos pelo wrapper de medição externa. Os demais contêm os tempos internos por rank coletados com `MPI_Wtime`.

O ensaio incompleto de P=6 não foi incluído, pois não pertence à série final e não possui o mesmo conjunto de repetições e medições externas.
