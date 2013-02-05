#$$first setup when new repo is created
git config --global user.name "Zhaoxia Meng"
git config --global user.email "Zhaoxia.Meng@cern.ch"
#no need to do above if aready setup

touch README.txt
git add README.txt
git commit -m "README please"
git push -u origin master
