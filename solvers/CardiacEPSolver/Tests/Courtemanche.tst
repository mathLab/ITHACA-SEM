<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Courtemanche Cell model</description>
    <executable>PrePacing</executable>
    <parameters>Courtemanche.xml</parameters>
    <files>
        <file description="Session File">Courtemanche.xml</file>
    </files>
    <metrics>
        <!--<metric type="Regex" id="1">
            <regex>^(\d*\.?\d*)\s*([+-]?\d*\.?\d*)</regex>
            <matches>
                <match>
                    <field>0.5</field>
                    <field>-81.0043</field>
                </match>
                <match>
                    <field>1</field>
                    <field>-81.0085</field>
                </match>
                <match>
                    <field>10</field>
                    <field>-2.33809</field>
                </match>
            </matches>
        </metric>-->
        <metric type="Regex" id="2">
            <regex>^#\s([\w]*)\s*([+-]?\d*\.?\d*)</regex>
            <matches>
                <match>
                    <field>u</field>
                    <field>-2.29798</field>
                </match>
                <match>
                    <field>m</field>
                    <field>0.99334</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>




