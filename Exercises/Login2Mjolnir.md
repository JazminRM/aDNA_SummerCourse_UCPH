# Connecting to the server

We will be working on an Amazon Web Server (AWS) that has all the software and data necessary for the course. 

The server consists of 64 CPUs (~2 per person) and 128Gb of RAM.

#### Participants usernames

We have created an account for each participant, your user name is your KU-ID and your password is your KU password.

## Mac, Linux, Ubuntu

1. Download the 'Key' file (```apgc-2021-key.pem.txt```) from the Dropbox link above, and save it in a local directory in your computer. We recommend having a dedicated directory where you can work during the course, since we will be uploading and downloading files from the server. 

2. Open the Terminal application (in a Mac computer you can find this on **Applications**, **Other**, **Terminal**)

3. Go to the directory where you have the Key file using ```cd```: 

```{bash, eval = FALSE}
cd /path/to/dir/with/key/file/
```

4. Change the permissions on the 'Key' file: 

```{bash, eval = FALSE}
chmod 400 apgc-2021-key.pem.txt
```

5. Login to the server by typing the following on your terminal (remember to change ```your_user_name``` for the username that was assigned to you):

```{bash, eval = FALSE}
ssh -i apgc-2021-key.pem.txt your_user_name@3.249.84.19
```

6. Finally, we will be downloading files from the server during the course, so check if you can download a test file from the server to your local computer with  ```scp ``` (remember to change ```your_user_name``` for the username that was assigned to you): 

```{bash, eval = FALSE}
scp -i apgc-2021-key.pem.txt your_user_name@3.249.84.19:/home/ec2-user/Data/testFile.txt .
```

And try uploading it to your home directory (this time there are two places where you should change ```your_user_name``` for your own username): 

```{bash, eval = FALSE}
scp -i apgc-2021-key.pem.txt testFile.txt your_user_name@3.249.84.19:/home/your_user_name/
```

Don't forget to exit by typing ```exit``` once you are done. 

Try login in and out a few times to make sure everything is ready for the course. 

## Windows

Windows users should download **PuTTY** to connect to their virtual server, and **WinSCP** to transfer files between their local computer and the server. 

1. Download the 'Key' file (```apgc-2021-key.pem.txt```) from the Dropbox link above, and save it in a local directory in your computer. We recommend having a dedicated directory where you can work during the course, since we will be uploading and downloading files from the server. 

2. Download and install PuTTY: https://www.chiark.greenend.org.uk/~sgtatham/putty/

3. Follow the instructions in sections **Convert your private key using PuTTYgen**, **Connect to your Linux instance** and **Transfer files to your Linux instance using WinSCP** of this AWS manual: https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/putty.html 

In step 2 of section 3 (Connect to your Linux instance) use the following information:

**my-instance-user-name**: the ```your_user_name``` that was assigned to you

**my-instance-public-dns-name**: 3.249.84.19

4. Finally, we will be downloading files from the server during the course, so once you installed and configured **WinSCP**, check if you can download this file to your computer:

```
/home/ec2-user/Data/testFile.txt 
```

And try uploading it again to your home directory:

```
/home/your_user_name
```
