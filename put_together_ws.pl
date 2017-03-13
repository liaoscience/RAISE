#!/usr/bin/perl
#根据路径把文件合并在一起,就是把一个文件夹里的文件合在一起.适合于基因表达谱之类的问题.
#如何读取文件夹, 如何push, hash数组. print 可以直接在一行后面加内容
#由王森创作的脚本，20150212
#

my $dir=$ARGV[0];
my @file;
my %hash;
opendir DIR,"$dir" or die $!; #打开路径
for my $file(readdir DIR) { #读取文件夹下的文件
        next unless $file=~/^\w+/;
        open IN,"$dir/$file" or die $!; #打开文件
        push @file,$file; #文件名放入数组
        while (<IN>) {
                chomp;
                my @array=split;
                #next unless ($array[1]=~/^\d+/);
                $hash{$array[0]}{$file}=$array[1]; #哪个是键，值 哈希数组。
        }
        close IN;
}
closedir DIR;
for my $i(keys %hash) {
        print "$i";
        for my $f(@file) {
                $c=($hash{$i}{$f})?$hash{$i}{$f}:0;
                print "\t$c";
        }
        print "\n";
}


