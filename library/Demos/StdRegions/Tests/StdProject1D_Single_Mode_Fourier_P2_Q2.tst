<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject1D Segment single mode Fourier basis P=2 Q=2</description>
    <executable>StdProject</executable>
    <parameters>-s segment -b FourierSingleMode -o 2 -p 2</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">2.77556e-17</value>
        </metric>
    </metrics>
</test>


