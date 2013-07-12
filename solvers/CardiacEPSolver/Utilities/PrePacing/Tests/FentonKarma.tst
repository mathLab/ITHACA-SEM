<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fenton-Karma cell model</description>
    <executable>PrePacing</executable>
    <parameters>FentonKarma.xml</parameters>
    <files>
        <file description="Session File">FentonKarma.xml</file>
    </files>
    <metrics>
        <metric type="Regex" id="1">
            <regex>
                ^#\s([\w]*)\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)
            </regex>
            <matches>
                <match>
                    <field>u</field>
                    <field tolerance="1e-06">1.18454</field>
                </match>
                <match>
                    <field>v</field>
                    <field tolerance="1e-06">1.6656e-13</field>
                </match>
                <match>
                    <field>w</field>
                    <field tolerance="1e-06">0.89349</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>




