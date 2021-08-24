<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff2D Quadrilateral Fourier basis P=6 Q=8</description>
    <executable>StdProject</executable>
    <parameters>-s quadrilateral -b Fourier Fourier -o 6 6 -p 8 8 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.54556e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">2.84217e-14</value>
        </metric>
    </metrics>
</test>


