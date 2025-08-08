## VPN

For the practical exercises we will be using our local computing cluster (**mjolnir**). To access **mjolnir**, you'll need to first login into the University of Copenhagen's VPN. 

If you don't have the VPN installed, you can download it and find the instructions on how to set it up for different OS here:

https://kunet.ku.dk/employee-guide/Pages/IT/Remote-access.aspx

## Connecting to mjolnir

Below you'll find the instruction on how to login to mjolnir for different OS. 

### Mac, Linux, Ubuntu

1. Login to the VPN (remember you'll need to accept the MFA request in your phone). 

2. Open the Terminal application (in a Mac computer you can find this on ```/Applications/Other/Terminal```).

3. Login to mjolnir by typing the following in your terminal (make sure to change ```ku_username``` for your own KU ID):

```{bash, eval = FALSE}
ssh ku_username@mjolnirgate.unicph.domain
# you'll get a prompt asking for your password, type it and hit <enter>
```
Don't forget to exit by typing ```exit``` once you are done running the exercises. 

#### Transfering files between mjolnir and your local computer

If you need to download/upload files from/to **mjolnir** you can do it through your terminal.

Download a file called "testFile.txt" to your current directory in your local computer:
```{bash, eval = FALSE}
# First let's check how the looks like: 
cat /projects/course_1/people/clx746/Data/testFile.txt 
```

Open a new terminal window and go into the directory where you want to download the file (in the example, I am downloading it in my Desktop): 
```{bash, eval = FALSE}
# go to the local directory: 
cd /Users/Jazmin/Desktop

# Then use scp to download the file (remember to change ku_username for your username/ku ID).  
scp ku_username@mjolnirgate.unicph.domain:/projects/course_1/people/clx746/Data/testFile.txt .
```

Open the file in your local computer with, for example, textedit and edit the file and save it. 

Now you can upload it back to the server like this:
```{bash, eval = FALSE}
# Remember to change ku_username for your username (ku ID).  
scp testFile.txt ku_username@mjolnirgate.unicph.domain:/projects/course_1/people/ku_username/
```

Finally, check if the file you uploaded changed: 
```{bash, eval = FALSE}
# Remember to change ku_username for your username (ku ID).  
cat /projects/course_1/people/ku_username/testFile.txt 
```

<p>&nbsp;</p>

### Windows

Windows users need to download **MobaXterm** (https://mobaxterm.mobatek.net/) to connect to the server, and **WinSCP** (https://winscp.net/eng/download.php) to transfer files between their local computer and the server. Download these two programs and follow the instructions for setting it up. 

Once you have WinSCP installed you'll be able to see and navigate within your local computer directories as well as the directories in mjolnir. Using WinSCP visual interface, go to the following directory and download the "testFile.txt". 

```
/projects/course_1/people/clx746/Data/testFile.txt 
```

Open the file in your local computer with, for example, WordPad or notePad and edit the file and save it. 

Using WinSCP visual interface again, try uploading the "testFile.txt" to your local directory: 

```
/projects/course_1/people/ku_username/
```

Try login in and out a few times to make sure everything works 







