/* Contains the functions used for the socket communication.

Contains both the functions that transmit data to the socket and read the data
back out again once finished, and the function which opens the socket initially.

Functions:
   error: Prints an error message and then exits.
   open_socket_: Opens a socket with the required host server, socket type and
      port number.
   write_buffer_: Writes a string to the socket.
   read_buffer_: Reads data from the socket.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h> 

void error(const char *msg)
// Prints an error message and then exits.

{   perror(msg);  exit(0);   }

void open_socket_(int *psockfd, int* inet, int* port, char* host)  
/* Opens a socket.

Note that fortran passes an extra argument for the string length, but this is 
ignored here for C compatibility.

Args:
   psockfd: The id of the socket that will be created.
   inet: An integer that determines whether the socket will be an inet or unix
      domain socket. Gives unix if 0, inet otherwise.
   port: The port number for the socket to be created. Low numbers are often
      reserved for important channels, so use of numbers of 4 or more digits is 
      recommended.
   host: The name of the host server.
*/

{
   int sockfd, portno, n;
   struct hostent *server;

   fprintf(stderr, "Connection requested %s, %d, %d\n", host, *port, *inet);
   struct sockaddr * psock; int ssock;
   if (*inet>0)
   {  
      struct sockaddr_in serv_addr;      psock=(struct sockaddr *)&serv_addr;     ssock=sizeof(serv_addr);
      sockfd = socket(AF_INET, SOCK_STREAM, 0);
      if (sockfd < 0)  error("ERROR opening socket");
   
      server = gethostbyname(host);
      if (server == NULL)
      {
         fprintf(stderr, "ERROR, no such host %s \n", host);
         exit(0);
      }

      bzero((char *) &serv_addr, sizeof(serv_addr));
      serv_addr.sin_family = AF_INET;
      bcopy((char *)server->h_addr, (char *)&serv_addr.sin_addr.s_addr, server->h_length);
      serv_addr.sin_port = htons(*port);
   }
   else
   {
      struct sockaddr_un serv_addr;      psock=(struct sockaddr *)&serv_addr;     ssock=sizeof(serv_addr);
      sockfd = socket(AF_UNIX, SOCK_STREAM, 0);
      bzero((char *) &serv_addr, sizeof(serv_addr));
      serv_addr.sun_family = AF_UNIX;
      strcpy(serv_addr.sun_path, "/tmp/wrappi_");
      strcpy(serv_addr.sun_path+12, host);
   }
   
   if (connect(sockfd, psock, ssock) < 0) error("ERROR connecting");

   *psockfd=sockfd;
}

void writebuffer_(int *psockfd, char *data, int* plen)
/* Writes to a socket.

Args:
   psockfd: The id of the socket that will be written to.
   data: The data to be written to the socket.
   plen: The length of the data in bytes.
*/

{
   int n;   
   int sockfd=*psockfd;
   int len=*plen;

   n = write(sockfd,data,len);
   if (n < 0) error("ERROR writing to socket");
}


void readbuffer_(int *psockfd, char *data, int* plen)
/* Reads from a socket.

Args:
   psockfd: The id of the socket that will be read from.
   data: The storage array for data read from the socket.
   plen: The length of the data in bytes.
*/

{
   int n, nr;
   int sockfd=*psockfd;
   int len=*plen;

   n = nr = read(sockfd,data,len);
   
   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n == 0) error("ERROR reading from socket");
}


