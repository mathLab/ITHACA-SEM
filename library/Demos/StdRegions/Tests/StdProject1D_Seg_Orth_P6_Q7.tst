<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject1D Segment Orthonormal basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s Segment -b Ortho_A -o 6 -p 7</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.37715e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.16573e-15</value>
        </metric>
    </metrics>
</test>


