<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject2D Quadrilateral Fourier Single Mode basis P=2 Q=2</description>
    <executable>StdProject</executable>
    <parameters>-s quadrilateral -b FourierSingleMode FourierSingleMode -o 2 2 -p 2 2</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">3.59678e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">8.88178e-16</value>
        </metric>
    </metrics>
</test>


