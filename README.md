# Preditor de Estrutura terciária


### pré-requisitos:

* Biblioteca Rosetta
* PyRosetta
* Sistema Operacional Linux
* gcc-9
* Credenciais do Labicom no Graylab:
* * user: levinthal
* * senha: paradox

****
**Ambiente Virtual**

É fortemente recomentado que voce utilize um ambiente virtual para instalação das dependências desse projeto, para isso instale o genrenciador de versões do python pyenv [seguindo este tutorial](https://k0nze.dev/posts/install-pyenv-venv-vscode/)

A versão recomentada do python para rodar este projeto é o python 3.10.12

Com essa versão selecionada a partir do pyenv, crie um ambiente virtual com o comando:

``python3 -m venv protein_prediction``

Recomenda-se que os próximos passos sejam feitos a partir do ambiente virtual.
___

**gcc-9.5**

Para compilar os arquivos binários do Rosetta é necessário que o seu compilador gcc e g++ estejam na versão 9, até o momento da escrita dessa documentação o Rosetta não suportava o gcc-10 ou gcc-11. Para alterar a versão do seu gcc, digite em seu terminal:

```sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 10```
```sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 10```

```sudo update-alternatives --set gcc /usr/bin/gcc-9```
```sudo update-alternatives --set g++ /usr/bin/g++-9```
___

**Instalar Rosetta** 

1. O primeiro passo para a instalação do Rosetta é acessar o link [para download](https://www.rosettacommons.org/software/academic);


2. Na tela de Downloads selecione a versão mais recente do Software;


3. Na opção de Download escolha a opção de source code;

Obs: O download pode demorar muito tempo a depender da velocidade de sua rede.

4. Após o Download ser concluído navegue até a pasta que o arquivo foi baixado e digite o comando:

``tar - xvzf <nome do arquivo baixado>.tgz``

5. Após descompactar o arquivo baixado, renomeio-o para Rosetta

6. Será necessária a instalação de alguns pacotes para instalar o Rosetta, para isso digite em seu terminal:

``sudo apt install zlib1g-dev scons build-essential``

7. Após o término da instalação dos poscotes necesários iremos de fato instalar o Rosetta, para isso, abra no terminal o diretório em que a pasta do Rosetta foi extraída e digite:

``cd Rosetta/main/source``

``./scons.py mode=release bin``

**Instalar PyRosetta** 

1. O primeiro passo para a instalação do PyRosetta é acessar o link [para download](https://graylab.jhu.edu/download/PyRosetta4/archive/release/) (use as credenciais se necessário);


2. Na tela de Downloads selecione a versã "PyRosetta4.Release.python310.ubuntu.cxx11thread.serialization/";


3. Na opção de Download escolha a opção lasted.html;

Obs: O download pode demorar muito tempo a depender da velocidade de sua rede.

4. Após o Download ser concluído navegue até a pasta que o arquivo foi baixado e digite o comando:

``tar - vjxf <nome do arquivo baixado>.tar.bz2``

5. Após descompactar o arquivo baixado, renomeio-o para PyRosetta

6. Pelo terminal, navegue até o diretório que a pasta foi descompactada e digite:

``python3 PyRosetta/setup/setup.py install``
