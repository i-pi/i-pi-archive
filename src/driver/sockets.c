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
{
   perror(msg);
   exit(0);
}

void open_socket_(int *psockfd)
{
   int sockfd, portno, n;
   struct hostent *server;

   struct sockaddr_in serv_addr;
   sockfd = socket(AF_INET, SOCK_STREAM, 0);
   if (sockfd < 0)
      error("ERROR opening socket");
   server = gethostbyname("spartaco.dyndns.org");
   if (server == NULL) 
   {
      fprintf(stderr, "ERROR, no such host\n");
      exit(0);
   }

   bzero((char *) &serv_addr, sizeof(serv_addr));
   serv_addr.sin_family = AF_INET;
   bcopy((char *)server->h_addr, 
      (char *)&serv_addr.sin_addr.s_addr,
      server->h_length);
   serv_addr.sin_port = htons(31415);

/*
   struct sockaddr_un serv_addr;
   sockfd = socket(AF_UNIX, SOCK_STREAM, 0);
   bzero((char *) &serv_addr, sizeof(serv_addr));
   serv_addr.sun_family = AF_UNIX;
   strcpy(serv_addr.sun_path, "/tmp/wrappi_localhost");
*/
   if (connect(sockfd,(struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) 
      error("ERROR connecting");

   *psockfd=sockfd;

}

void writebuffer_(int *psockfd, char *data, int* plen)
{
   int n;   
   int sockfd=*psockfd;
   int len=*plen;

   n = write(sockfd,data,len);
   if (n < 0) 
      error("ERROR writing to socket");
}


void readbuffer_(int *psockfd, char *data, int* plen)
{
   int n, nr;
   int sockfd=*psockfd;
   int len=*plen;

   n = nr = read(sockfd,data,len);
   
   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n == 0)
      error("ERROR reading from socket");
}


