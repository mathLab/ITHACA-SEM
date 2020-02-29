<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>BaryInterp2D Triangle Modified basis P=6 Q=7</description>
    <executable>BaryInterp</executable>
    <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7</parameters>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-12">0.0124512</value>
        </metric>
        <metric type="L2" id="2">
            <value tolerance="1e-12">0.00715705</value>
        </metric>
    </metrics>
</test>


