### Instructions for setting up the course on CoCalc

1. Login to CoCalc.
2. If you haven't already done so, setup your [ssh keys](https://cocalc.com/settings/ssh-keys?session=default).
3. Create a new project.
4. Enter that project.
5. Create a course file (our gitignore expects the prefix to be GCB535-anything). We use GCB535-YYYY.
6. Setup deploy keys:
    * Click on project settings. Look for SSH keys.
    * ssh to the username@host specified in the section. This should work based on the account ssh key setup that you did in 2.
    * (Aside: You'll need both rsa and dsa deploy keys if you want to pushes changes to both the main repo and the answer keys repo. The following instructions setup both.)
    * Once you're in via ssh, run `ssh-keygen`. This will generate an rsa key. The password you use here will be required to push/pull changes from your github repo.
    * Run `cat ~/.ssh/id_rsa.pub`. This will show the public key. On GitHub, navigate to your GCB535 repo, settings, deploy keys. Do "add deploy key" and paste the public key. Select "allow write access" if you intend to push changes back.
    * Run run `ssh-keygen /home/user/.ssh/id_rsa_anskey`. This will generate an rsa key for the answer key. The password you use here will be required to push/pull changes from your github repo.
    * Run `cat ~/.ssh/id_rsa_anskey.pub`. This will show the public key. On GitHub, navigate to your GCB535-Keys repo, settings, deploy keys. Do "add deploy key" and paste the public key. Select "allow write access" if you intend to push changes back.
7. Run git clone for the main GCB535 repo. (for us: `git clone git@github.com:greenelab/GCB535.git`)
8. Run `cd GCB535`
9. Run `mv * ../`
10. Run `mv .* ../` - you'll get some errors about . and .. directories, but can ignore these. The key is to move the git info up a level.
11. Cleanup: Run `cd ..` and then `rm -r GCB535`. This will get the exercises into the root folder, which makes the CoCalc class functionality happy.
12. If you followed the instructions above and your answer keys use the answer key rsa key, we're going to make keygit an alias for git using the DSA key. Run: `echo "alias keygit=\"GIT_SSH_COMMAND='ssh -i ~/.ssh/id_rsa_anskey' git\"" >> ~/.bashrc`
13. run `source .bashrc`
14. Get the answer keys onto SMC as well. Use the keygit alias to clone your answer key repo (for us: `keygit clone git@github.com:greenelab/GCB535-Keys.git`). NOTE: TO DO ANYTHING WITH GIT ON THE ANSWER KEYS, USE THE `keygit` ALIAS.

Congratulations! Your GitHub setup tasks are done!
